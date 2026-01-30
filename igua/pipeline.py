import abc
import dataclasses
import pathlib
import tempfile
import typing
from typing import Dict, Literal

import anndata
import rich.console
import rich.progress
import scipy.sparse
import numpy
import pandas

from .dataset.base import BaseDataset
from .dataset.defensefinder import DefenseFinderDataset
from .dataset.fasta_gff import FastaGFFDataset
from .mmseqs import MMSeqs, Database
from .clustering import HierarchicalClustering, LinearClustering


class _BaseSink(abc.ABC):

    def __init__(self):
        self.stats = []
        self.done = set()

    def add_record(self, name: str, sequence: str, **kwargs) -> bool:
        if name not in self.done:
            self.stats.append({"id": name, "length": len(sequence), **kwargs})
            self.done.add(name)
            return True
        else:
            return False

    def report_statistic(self) -> pandas.DataFrame:
        return pandas.DataFrame(self.stats)


class _FASTASink(_BaseSink):

    def __init__(self, file: typing.TextIO) -> None:
        super().__init__()
        self.file = file

    def add_record(self, name: str, sequence: str, **kwargs) -> None:
        if not super().add_record(name, sequence, **kwargs):
            return False
        self.file.write(">")
        self.file.write(name)
        self.file.write("\n")
        self.file.write(sequence)
        self.file.write("\n")
        return True


@dataclasses.dataclass
class ClusteringParameters:
    nuc1: Dict[str, object]
    nuc2: Dict[str, object]
    prot: Dict[str, object]
    clustering_method: Literal["average", "single", "complete", "weighted", "centroid", "median", "ward"]
    clustering_distance: float
    precision: Literal["half", "single", "double"]

    @classmethod
    def default(cls) -> "ClusteringParameters":
        """Create new default clustering parameters.
        """
        return cls(
            nuc1=dict(
                e_value=0.001,
                sequence_identity=0.85,
                coverage=1.0,
                cluster_mode=0,
                coverage_mode=1,
                spaced_kmer_mode=0,
            ),
            nuc2=dict(
                e_value=0.001,
                sequence_identity=0.6,
                coverage=0.5,
                cluster_mode=0,
                coverage_mode=0,
                spaced_kmer_mode=0,
            ),
            prot=dict(
                e_value=0.001,
                coverage=0.9,
                coverage_mode=1,
                sequence_identity=0.5,
            ),
            clustering_method="average",
            clustering_distance=0.8,
            precision="double",
        )


@dataclasses.dataclass
class ClusteringResult:
    gcfs: pandas.DataFrame
    compositions: typing.Optional[anndata.AnnData]
    proteins_faa: typing.Optional[pathlib.Path]


class ClusteringPipeline:
    """The IGUA multi-stage clustering pipeline.
    """

    def __init__(
        self,
        workdir: pathlib.Path,
        params: typing.Optional[ClusteringParameters] = None,
        *,
        prefix: str = "GCF",
        jobs: int = 1,
        mmseqs: typing.Optional[MMSeqs] = None,
        progress: typing.Optional[rich.progress.Progress] = None,
    ):
        self.jobs = jobs
        self.params = params or ClusteringParameters.default()
        self.workdir = pathlib.Path(workdir)
        self.prefix = prefix

        if mmseqs is None:
            self.mmseqs = MMSeqs(progress=progress, threads=self.jobs, tempdir=self.workdir)
        else:
            self.mmseqs = mmseqs

        self.progress = progress
        if progress is None:
            self.console = rich.console.Console(quiet=True)
        else:
            self.console = progress.console

        if self.params.clustering_method == "linclust":
            self.clustering = LinearClustering(
                distance=self.params.clustering_distance,
            )
        else:
            self.clustering = HierarchicalClustering(
                method=self.params.clustering_method,
                distance=self.params.clustering_distance,
                precision=self.params.precision,
                jobs=self.jobs,
            )

    # ---

    def extract_clusters_to_file(
        self,
        dataset: BaseDataset,
        output: pathlib.Path,
    ):
        # record processed sequences to flag duplicates
        done = set()
        n_duplicates = 0
        # extract raw sequences
        self.console.print(f"[bold blue]{'Loading':>12}[/] input clusters")
        with output.open("w") as dst:
            sink = _FASTASink(dst)
            for cluster in dataset.extract_clusters(self.progress):
                if cluster.id not in done:
                    done.add(cluster.id)
                    sink.add_record(cluster.id, cluster.sequence, filename=cluster.source)
                else:
                    n_duplicates += 1
            input_sequences = (
                sink.report_statistic()
                    .rename(columns={"id": "cluster_id", "length": "cluster_length"})
                    .set_index("cluster_id")
            )
        if n_duplicates > 0:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] {n_duplicates} duplicate "
                "clusters found in inputs"
            )

        self.console.print(
            f"[bold green]{'Loaded':>12}[/] {len(input_sequences)} "
            "clusters to process"
        )
        return input_sequences

    def extract_proteins_to_file(
        self,
        dataset: BaseDataset,
        clusters: typing.Container[str],
        output: pathlib.Path,
    ) -> pandas.DataFrame:
        self.console.print(
            f"[bold blue]{'Extracting':>12}[/] protein sequences from "
            "representative clusters"
        )
        with output.open("w") as dst:
            sink = _FASTASink(dst)
            for protein in dataset.extract_proteins(self.progress, clusters):
                sink.add_record(protein.id, protein.sequence, cluster_id=protein.cluster_id)
            protein_sizes = (
                sink.report_statistic()
                    .rename(columns={"id": "protein_id", "length": "protein_length"})
                    .set_index("protein_id")
            )
        return protein_sizes

    # ---

    def make_compositions(
        self,
        protein_clusters: pandas.DataFrame,
        representatives: typing.Dict[str, int],
        protein_representatives: typing.Dict[str, int],
    ) -> anndata.AnnData:
        # index the proteins by protein ID
        lengths = protein_clusters.set_index("protein_id")["protein_length"]

        # create the empty compositional matrix
        compositions = scipy.sparse.dok_matrix(
            (len(representatives), len(protein_representatives)), 
            dtype=numpy.int32
        )

        # Build compositional matrix using protein lengths as weights
        task = self.progress.add_task(
            description=f"[bold blue]{'Working':>9}[/]", 
            total=len(protein_clusters)
        )
        for row in self.progress.track(protein_clusters.itertuples(), task_id=task):
            cluster_index = representatives[row.cluster_id]
            prot_index = protein_representatives[row.protein_representative]
            compositions[cluster_index, prot_index] += lengths.loc[row.protein_representative]
        self.progress.remove_task(task)

        # Build the annotated data matrix
        sorted_representatives = sorted(representatives, key=representatives.__getitem__)
        sorted_protein_representatives = sorted(
            protein_representatives, key=protein_representatives.__getitem__
        )
        return anndata.AnnData(
            X=compositions.tocsr(),
            obs=pandas.DataFrame(
                index=pandas.Index(sorted_representatives, name="cluster_id")
            ),
            var=pandas.DataFrame(
                index=pandas.Index(sorted_protein_representatives, name="protein_id"),
                data=dict(size=lengths.loc[sorted_protein_representatives]),
            ),
        )

    # ---

    def run(
        self,
        dataset: BaseDataset,
        *,
        clustering: bool = True,
    ):
        # extract raw sequences
        clusters_fna = self.workdir / "clusters.fna"
        input_sequences = self.extract_clusters_to_file(dataset, output=clusters_fna)

        # create initial sequence database
        self.console.print(
            f"[bold blue]{'Starting':>12}[/] nucleotide deduplication step "
            "with [purple]MMSeqs2[/]"
        )
        db = Database.create(self.mmseqs, clusters_fna)
        step1 = db.cluster(self.workdir / "step1.db", **self.params.nuc1)
        gcfs1 = step1.to_dataframe(columns=["fragment_representative", "cluster_id"]).sort_values("cluster_id")  # type: ignore
        self.console.print(
            f"[bold green]{'Reduced':>12}[/] {len(gcfs1)} clusters to {gcfs1.fragment_representative.nunique()} complete representatives"
        )

        # cluster sequences
        self.console.print(
            f"[bold blue]{'Starting':>12}[/] nucleotide clustering step with [purple]MMSeqs2[/]"
        )
        repdb = step1.to_subdb(self.workdir / "step1.rep_seq.db")
        step2 = repdb.cluster(self.workdir / "step2.db", **self.params.nuc2)
        gcfs2 = step2.to_dataframe(columns=["nucleotide_representative", "fragment_representative"]).sort_values("fragment_representative")  # type: ignore
        self.console.print(
            f"[bold green]{'Reduced':>12}[/] {len(gcfs2)} clusters to {len(gcfs2.nucleotide_representative.unique())} nucleotide representatives"
        )

        # load representatives
        self.console.print(
            f"[bold blue]{'Extracting':>12}[/] representative clusters"
        )
        representatives = {
            x: i
            for i, x in enumerate(sorted(gcfs2["nucleotide_representative"].unique()))
        }
        self.console.print(
            f"[bold green]{'Found':>12}[/] {len(representatives)} nucleotide representative clusters"
        )

        if clustering and len(representatives) > 1:
            # extract proteins and record sizes
            proteins_faa = self.workdir / "proteins.faa"
            protein_sizes = self.extract_proteins_to_file(dataset, representatives, output=proteins_faa)
            if not proteins_faa.exists() or proteins_faa.stat().st_size == 0:
                self.console.print(
                    f"[bold yellow]{'Warning':>12}[/] No proteins extracted from input dataset"
                )

            # cluster proteins
            prot_db = Database.create(self.mmseqs, proteins_faa)
            prot_result = prot_db.cluster(self.workdir / "step3.db", **self.params.prot)
            prot_clusters = prot_result.to_dataframe(
                columns=["protein_representative", "protein_id"]
            )

            # record protein lengths and source cluster
            prot_clusters = pandas.merge(
                prot_clusters, 
                protein_sizes[["cluster_id", "protein_length"]], 
                on="protein_id"
            )
       
            # extract protein representatives
            protein_representatives = {
                x: i
                for i, x in enumerate(
                    sorted(prot_clusters["protein_representative"].unique())
                )
            }
            self.console.print(
                f"[bold green]{'Found':>12}[/] {len(protein_representatives)} "
                f"protein representatives for {len(prot_clusters)} proteins"
            )

            # build weighted compositional array
            self.console.print(
                f"[bold blue]{'Building':>12}[/] weighted compositional array"
            )
            compositions = self.make_compositions(
                prot_clusters,
                representatives,
                protein_representatives,
            )

            # run clustering based on protein array membership
            self.console.print(
                f"[bold blue]{'Clustering':>12}[/] gene clusters using "
                f"{self.params.clustering_method} linkage"
            )
            flat = self.clustering.cluster(compositions.X)

            # build GCFs based on flat clustering
            gcfs3 = pandas.DataFrame(
                {
                    "gcf_id": [f"{self.prefix}{i:07}" for i in flat],
                    "nucleotide_representative": compositions.obs_names,
                }
            )
        
        else:
            compositions = None
            proteins_faa = None
            sorted_representatives = sorted(
                representatives, key=representatives.__getitem__
            )
            gcfs3 = pandas.DataFrame(
                {
                    "gcf_id": [
                        f"{self.prefix}{i+1:07}"
                        for i in range(len(sorted_representatives))
                    ],
                    "nucleotide_representative": sorted_representatives,
                }
            )

        self.console.print(
            f"[bold green]{'Built':>12}[/] {len(gcfs3.gcf_id.unique())} "
            f"GCFs from {len(input_sequences)} initial clusters"
        )

        # build GCFs based on flat clustering
        gcf3_representatives = (
            pandas.merge(
                gcfs3,
                input_sequences["cluster_length"],
                left_on="nucleotide_representative",
                right_index=True,
            )
            .sort_values("cluster_length")
            .drop_duplicates("gcf_id", keep="last")
            .set_index("gcf_id")
        )
        gcfs3 = pandas.merge(
            gcfs3,
            gcf3_representatives["nucleotide_representative"].rename(
                "gcf_representative"
            ),
            left_on="gcf_id",
            right_index=True,
        )

        # build final GCF table
        gcfs = pandas.merge(
            pandas.merge(
                pandas.merge(gcfs1, gcfs2, on="fragment_representative"),
                gcfs3,
                on="nucleotide_representative",
            ),
            input_sequences,
            left_on="cluster_id",
            right_index=True,
        )

        # sort results and drop unused columns
        gcfs.sort_values(["gcf_id", "cluster_length"], inplace=True)
        gcfs = gcfs[
            [
                "cluster_id",
                "cluster_length",
                "gcf_id",
                "gcf_representative",
                "nucleotide_representative",
                "fragment_representative",
                "filename",
            ]
        ]
        
        return ClusteringResult(
            gcfs=gcfs, 
            compositions=compositions,
            proteins_faa=proteins_faa,
        )
