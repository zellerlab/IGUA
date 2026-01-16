import abc
import gzip
import io
import pathlib
import typing
from io import StringIO

import Bio.Seq
import gb_io
import pandas as pd
import polars as pl
import rich.progress

from .cluster_extractor import (
    ClusterDataAdapter,
    ClusterMetadataCache,
    GeneClusterExtractor,
    GenomeContext,
)
from .mmseqs import Database, MMSeqs


_GZIP_MAGIC = b"\x1f\x8b"


class BaseDataset(abc.ABC):
    """Base class for dataset extraction.
    This class defines the basic structure and methods for extracting nucleotide and
    protein sequences from various file formats. It serves as a base class for specific
    dataset classes like GenBankDataset and GFFDataset.
    """

    @abc.abstractmethod
    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
    ) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        pass

    def write_fasta(self, file: typing.TextIO, name: str, sequence: str) -> None:
        file.write(">{}\n".format(name))
        file.write(sequence)
        file.write("\n")
        return None

    def translate_orf(
        self, sequence: typing.Union[str, bytes], translation_table: int = 11
    ) -> str:
        return str(Bio.Seq.Seq(sequence).translate(translation_table))

    def create_sequence_database(
        self,
        mmseqs: MMSeqs,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output_db_path: pathlib.Path,
    ) -> Database:
        """Default implementation creates a temporary file then a database"""
        tmp_fasta = output_db_path.with_suffix(".fna")
        self.extract_sequences(progress, inputs, tmp_fasta)
        return Database.create(mmseqs, tmp_fasta)

    def create_protein_database(
        self,
        mmseqs: MMSeqs,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        representatives: typing.Container[str],
        output_db_path: pathlib.Path,
    ) -> typing.Tuple[Database, typing.Dict[str, int]]:
        """Default implementation creates a temporary file then a database"""
        tmp_fasta = output_db_path.with_suffix(".faa")
        protein_sizes = self.extract_proteins(
            progress, inputs, tmp_fasta, representatives
        )
        return Database.create(mmseqs, tmp_fasta), protein_sizes


class GenBankDataset(BaseDataset):
    """GenBank dataset class.
    This class is used to extract nucleotide and protein sequences from GenBank files.
    It inherits from the BaseDataset class and implements the extract_sequences
    and extract_proteins methods.
    """

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
    ) -> pd.DataFrame:
        """Extracts nucleotide sequences from GenBank files.
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted sequences.
        Returns:
            pd.DataFrame: DataFrame containing the extracted sequences.
        """
        data = []
        done = set()
        n_duplicate = 0
        with open(output, "w") as dst:
            task1 = progress.add_task(f"[bold blue]{'Working':>9}[/]")
            for input_path in progress.track(inputs, task_id=task1):
                task2 = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
                with io.BufferedReader(progress.open(input_path, "rb", task_id=task2)) as reader:  # type: ignore
                    if reader.peek().startswith(_GZIP_MAGIC):
                        reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
                    for record in gb_io.iter(reader):
                        if record.name in done:
                            n_duplicate += 1
                        else:
                            self.write_fasta(
                                dst, record.name, record.sequence.decode("ascii")
                            )
                            data.append((record.name, len(record.sequence), input_path))
                            done.add(record.name)
                progress.remove_task(task2)
            progress.remove_task(task1)
        if n_duplicate > 0:
            progress.console.print(
                f"[bold yellow]{'Skipped':>12}[/] {n_duplicate} clusters with duplicate identifiers"
            )
        return pd.DataFrame(
            data=data, columns=["cluster_id", "cluster_length", "filename"]
        ).set_index("cluster_id")

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Extracts protein sequences from GenBank files.
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted protein sequences.
            representatives (typing.Container[str]): Set of representative cluster IDs.
        Returns:
            typing.Dict[str, int]: Dictionary containing protein IDs and their sizes.
        """
        protein_sizes = {}
        with output.open("w") as dst:
            for input_path in inputs:
                task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
                with io.BufferedReader(progress.open(input_path, "rb", task_id=task)) as reader:  # type: ignore
                    if reader.peek()[:2] == b"\x1f\x8b":
                        reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
                    for record in gb_io.iter(reader):
                        if record.name in representatives:
                            for i, feat in enumerate(
                                filter(lambda f: f.kind == "CDS", record.features)
                            ):
                                qualifier = next(
                                    (
                                        qualifier
                                        for qualifier in feat.qualifiers
                                        if qualifier.key == "translation"
                                    ),
                                    None,
                                )
                                if qualifier is None:
                                    rich.print(
                                        f"[bold yellow]{'Warning':>12}[/] no 'translation' qualifier found in CDS feature of {record.name!r}"
                                    )
                                    translation = self.translate_orf(
                                        record.sequence[
                                            feat.location.start : feat.location.end
                                        ]
                                    )
                                else:
                                    translation = qualifier.value.rstrip("*")
                                protein_id = "{}_{}".format(record.name, i)
                                if protein_id not in protein_sizes:
                                    self.write_fasta(dst, protein_id, translation)
                                    protein_sizes[protein_id] = len(translation)
                progress.remove_task(task)
        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(protein_sizes):,} proteins from {len(representatives):,} nucleotide representative"
        )
        return protein_sizes


class FastaGFFDataset(BaseDataset):
    """FastaGFF dataset class.

    Uses adapters to handle different cluster data formats.
    """

    def __init__(self, adapter: ClusterDataAdapter) -> None:
        """Initialize the FastaGFFDataset class.

        Args:
            adapter: Adapter for handling format-specific logic.
        """
        self.cluster_metadata: typing.Optional[typing.Union[pathlib.Path, str]] = None
        self.verbose: bool = False
        self.adapter = adapter
        self.gff_cache_dir: typing.Optional[pathlib.Path] = None
        self._metadata_cache_path: typing.Optional[pathlib.Path] = None

    def _create_genome_context(self, row: typing.Dict, genome_id: str) -> GenomeContext:
        """Create GenomeContext from a dataframe row."""
        return GenomeContext(
            genome_id=genome_id,
            systems_tsv=pathlib.Path(row["systems_tsv"]),
            genes_tsv=pathlib.Path(row["genes_tsv"]),
            gff_file=pathlib.Path(row["gff_file"]),
            genomic_fasta=pathlib.Path(row["genome_fasta_file"]),
            protein_fasta=pathlib.Path(row["protein_fasta_file"]),
            adapter=self.adapter,
        )

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
    ) -> pd.DataFrame:
        """Extracts nucleotide sequences from gene clusters and caches metadata."""
        extractor = GeneClusterExtractor(progress=progress, verbose=self.verbose)
        progress.console.print(
            f"[bold blue]{'Using':>12}[/] cluster metadata file: [magenta]{self.cluster_metadata}[/]"
        )

        df = pl.read_csv(self.cluster_metadata, separator="\t")
        self._metadata_cache_path = output.parent / f".{output.stem}_metadata.json"
        return self._extract_and_cache_metadata(progress, df, output, extractor)


    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Extracts protein sequences using cached metadata."""
        extractor = GeneClusterExtractor(progress=progress, verbose=self.verbose)

        if self._metadata_cache_path and self._metadata_cache_path.exists():
            progress.console.print(
                f"[bold blue]{'Using':>12}[/] cached metadata from sequence extraction"
            )
            return self._extract_proteins_from_cache(
                progress, output, extractor, representatives
            )

        progress.console.print(
            f"[bold yellow]{'Warning':>12}[/] No cached metadata found, reading cluster metadata again"
        )
        progress.console.print(
            f"[bold blue]{'Using':>12}[/] cluster metadata file: [magenta]{self.cluster_metadata}[/]"
        )

        try:
            df = pl.read_csv(self.cluster_metadata, separator="\t")
            return self._extract_proteins_direct(
                progress, df, output, extractor, representatives
            )
        except Exception as e:
            progress.console.print(
                f"[bold red]{'Error':>12}[/] reading cluster metadata: {e}"
            )
            return {}

    def _extract_and_cache_metadata(
        self,
        progress: rich.progress.Progress,
        df: pl.DataFrame,
        output: pathlib.Path,
        extractor: GeneClusterExtractor,
    ) -> pd.DataFrame:
        """Extract sequences and cache all metadata for protein extraction."""
        results = []

        def metadata_generator():
            """Generator that yields metadata for each genome."""
            genome_count = 0
            task = progress.add_task(
                f"[bold blue]{'Processing':>9}[/] gene clusters", total=len(df)
            )

            with open(output, "w") as dst:
                for row in df.iter_rows(named=True):
                    genome_id = row.get("genome_id") or f"genome_{genome_count:07}"
                    genome_count += 1
                    progress.update(
                        task,
                        description=f"[bold blue]{'Processing':>9}[/] strain: [bold cyan]{genome_id}",
                    )

                    context = self._create_genome_context(row, genome_id)

                    if not context.is_valid():
                        progress.console.print(
                            f"[bold yellow]{'Missing':>12}[/] files for {genome_id}"
                        )
                        progress.update(task, advance=1)
                        continue

                    try:
                        coordinates = extractor.extract_systems(context, dst)

                        for coord in coordinates:
                            if coord.valid:
                                results.append(
                                    (
                                        coord.sys_id,
                                        coord.end_coord - coord.start_coord + 1,
                                        coord.fasta_file,
                                    )
                                )

                        if coordinates:
                            yield {
                                "genome_id": genome_id,
                                "protein_fasta": str(context.protein_fasta),
                                "coordinates": [c.to_dict() for c in coordinates],
                            }

                    except Exception as e:
                        progress.console.print(
                            f"[bold red]{'Error':>12}[/] processing {genome_id}: {e}"
                        )

                    progress.update(task, advance=1)

            progress.remove_task(task)

        cache = ClusterMetadataCache(self._metadata_cache_path)
        cache.save_streaming(metadata_generator(), len(df))

        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(results):,} gene clusters from {len(df):,} strains/genomes"
        )
        progress.console.print(
            f"[bold blue]{'Cached':>12}[/] metadata to [magenta]{self._metadata_cache_path.name}[/]"
        )

        return pd.DataFrame(
            data=results, columns=["cluster_id", "cluster_length", "filename"]
        ).set_index("cluster_id")

    def _extract_proteins_from_cache(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
        extractor: GeneClusterExtractor,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Extract proteins using cached metadata with streaming."""
        cache = ClusterMetadataCache(self._metadata_cache_path)

        genome_count = sum(1 for _ in cache.iter_genomes())

        protein_sizes = {}

        with open(output, "w") as dst:
            task = progress.add_task(
                f"[bold blue]{'Processing':>9}[/] protein sequences", total=genome_count
            )

            for genome_metadata in cache.iter_genomes():
                genome_id = genome_metadata["genome_id"]
                progress.update(
                    task,
                    description=f"[bold blue]{'Processing':>9}[/] strain: [bold cyan]{genome_id}",
                )

                try:
                    proteins = extractor.extract_proteins_from_metadata(
                        genome_metadata, dst, representatives
                    )
                    protein_sizes.update(proteins)
                except Exception as e:
                    progress.console.print(
                        f"[bold red]{'Error':>12}[/] processing {genome_id}: {e}"
                    )

                progress.update(task, advance=1)

            progress.remove_task(task)

        self._log_protein_summary(progress, protein_sizes, representatives)
        return protein_sizes

    def _extract_proteins_direct(
        self,
        progress: rich.progress.Progress,
        df: pl.DataFrame,
        output: pathlib.Path,
        extractor: GeneClusterExtractor,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Direct protein extraction without cache (fallback)."""
        protein_sizes = {}

        with open(output, "w") as dst:
            task = progress.add_task(
                f"[bold blue]{'Processing':>9}[/] protein sequences", total=len(df)
            )

            for row in df.iter_rows(named=True):
                genome_id = row.get("genome_id") or f"genome_{0:07}"
                progress.update(
                    task,
                    description=f"[bold blue]{'Processing':>9}[/] strain: [bold cyan]{genome_id}",
                )

                context = self._create_genome_context(row, genome_id)

                if not context.is_valid():
                    progress.console.print(
                        f"[bold yellow]{'Missing':>12}[/] files for {genome_id}"
                    )
                    progress.update(task, advance=1)
                    continue

                try:
                    coordinates = []
                    temp_buffer = StringIO()
                    coords = extractor.extract_systems(context, temp_buffer)

                    genome_metadata = {
                        "genome_id": genome_id,
                        "protein_fasta": str(context.protein_fasta),
                        "coordinates": [c.to_dict() for c in coords],
                    }

                    proteins = extractor.extract_proteins_from_metadata(
                        genome_metadata, dst, representatives
                    )
                    protein_sizes.update(proteins)

                except Exception as e:
                    progress.console.print(
                        f"[bold red]{'Error':>12}[/] processing {genome_id}: {e}"
                    )

                progress.update(task, advance=1)

            progress.remove_task(task)

        self._log_protein_summary(progress, protein_sizes, representatives)
        return protein_sizes

    def _log_protein_summary(
        self,
        progress: rich.progress.Progress,
        results: typing.Dict,
        representatives: typing.Optional[typing.Container[str]],
    ):
        """Log summary for protein extraction."""
        rep_count = "all"
        if representatives:
            try:
                rep_count = str(len(representatives))
            except TypeError:
                rep_count = "specified"

        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(results):,} proteins from {rep_count} representative gene clusters"
        )
