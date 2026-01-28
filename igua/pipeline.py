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
from scipy.cluster.hierarchy import fcluster

from .dataset.base import BaseDataset
from .dataset.defensefinder import DefenseFinderDataset
from .dataset.fasta_gff import FastaGFFDataset
from .mmseqs import MMSeqs, Database
from .hca import manhattan, linkage


@dataclasses.dataclass
class ClusteringParameters:
    nuc1: Dict[str, object]  # TODO
    nuc2: Dict[str, object]  # TODO
    prot: Dict[str, object]  # TODO
    clustering_method: Literal["average", "single", "complete", "weighted", "centroid", "median", "ward"]
    clustering_distance: float


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
        params: ClusteringParameters,
        workdir: pathlib.Path,
        *,
        prefix: str = "GCF",
        jobs: int = 1,
        mmseqs: typing.Optional[MMSeqs] = None,
        progress: typing.Optional[rich.progress.Progress] = None,
        precision: Literal["half", "single", "double"] = "double",
    ):
        self.jobs = jobs
        self.params = params
        self.workdir = pathlib.Path(workdir)
        self.precision = precision
        self.prefix = prefix

        if mmseqs is None:
            self.mmseqs = MMSeqs(progress=progress, threads=jobs, tempdir=self.workdir)
        else:
            self.mmseqs = mmseqs

        self.progress = progress
        if progress is None:
            self.console = rich.console.Console(quiet=True)
        else:
            self.console = progress.console

    # ---

    def compute_distances(
        self,
        compositions: scipy.sparse.csr_matrix,
        jobs: typing.Optional[int],
        precision: str,
    ) -> numpy.ndarray:
        n = 0
        r = compositions.shape[0]
        # compute the number of amino acids per cluster
        clusters_aa = numpy.zeros(r, dtype=numpy.int32)
        clusters_aa[:] = compositions.sum(axis=1).A1
        # make sure the sparse matrix has sorted indices (necessary for
        # the distance algorithm to work efficiently)
        if not compositions.has_sorted_indices:
            compositions.sort_indices()
        # compute manhattan distance on sparse matrix
        distance_vector = numpy.zeros(r * (r - 1) // 2, dtype=precision)
        manhattan(
            compositions.data,
            compositions.indices,
            compositions.indptr,
            distance_vector,
            threads=jobs,
        )
        # ponderate by sum of amino-acid distance
        for i in range(r - 1):
            l = r - (i + 1)
            distance_vector[n : n + l] /= (clusters_aa[i + 1 :] + clusters_aa[i]).clip(
                min=1
            )
            n += l
        # check distances are in [0, 1]
        return numpy.clip(distance_vector, 0.0, 1.0, out=distance_vector)

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
        # check mmseqs version
        try:
            v = self.mmseqs.version()
            self.console.print(f"[bold green]{'Using':>12}[/] MMseqs2 {v!r}")
        except RuntimeError:
            self.console.print(
                f"[bold red]{'Failed':>12}[/] to locate MMseqs2 binary"
            )
            return errno.ENOENT

        # extract raw sequences
        clusters_fna = self.workdir / "clusters.fna"
        self.console.print(f"[bold blue]{'Loading':>12}[/] input clusters")
        input_sequences = dataset.extract_sequences(self.progress, clusters_fna)
        self.console.print(
            f"[bold green]{'Loaded':>12}[/] {len(input_sequences)} clusters to process"
        )

        # create initial sequence database
        self.console.print(
            f"[bold blue]{'Starting':>12}[/] nucleotide deduplication step with [purple]mmseqs[/]"
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
            self.console.print(
                f"[bold blue]{'Extracting':>12}[/] protein sequences from representative clusters"
            )

            protein_sizes = dataset.extract_proteins(
                self.progress, proteins_faa, representatives
            )

            if not proteins_faa.exists() or proteins_faa.stat().st_size == 0:
                self.console.print(
                    f"[bold yellow]{'Warning':>12}[/] No proteins extracted from defense systems"
                )
                self.console.print(
                    f"[bold yellow]{'Skipping':>12}[/] protein clustering due to empty protein file"
                )
                args.clustering = False

            # cluster proteins
            prot_db = Database.create(self.mmseqs, proteins_faa)
            prot_result = prot_db.cluster(self.workdir / "step3.db", **self.params.prot)
            prot_clusters = prot_result.to_dataframe(
                columns=["protein_representative", "protein_id"]
            )

            # record protein lengths and source cluster
            prot_clusters = pandas.merge(
                prot_clusters, 
                protein_sizes[["cluster_id", "protein_id", "protein_length"]], 
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

            # compute and ponderate distances
            self.console.print(
                f"[bold blue]{'Computing':>12}[/] pairwise distance based on "
                "protein composition"
            )
            distance_vector = self.compute_distances(
                compositions.X, 
                self.jobs, 
                self.precision
            )

            # run hierarchical clustering
            self.console.print(
                f"[bold blue]{'Clustering':>12}[/] gene clusters using "
                f"{self.params.clustering_method} linkage"
            )
            Z = linkage(distance_vector, method=self.params.clustering_method)
            flat = fcluster(Z, criterion="distance", t=self.params.clustering_distance)

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
