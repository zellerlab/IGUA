import pathlib
import typing
from io import StringIO


import pandas as pd
import polars as pl
import rich.progress

from .base_dataset import BaseDataset
from .cluster_extractor import (
    ClusterDataAdapter,
    ClusterMetadataCache,
    GeneClusterExtractor,
    GenomeContext,
)


class FastaGFFDataset(BaseDataset):
    """FastaGFF dataset class."""

    def __init__(
        self,
        inputs: typing.List[pathlib.Path],
        adapter: ClusterDataAdapter,
    ) -> None:
        """Initialize the FastaGFFDataset class.

        Args:
            inputs: List of input file paths (typically metadata TSV).
            adapter: Adapter for handling format-specific logic.
        """
        super().__init__(inputs)
        self.adapter = adapter
        self.verbose: bool = False
        self.gff_cache_dir: typing.Optional[pathlib.Path] = None
        self._metadata_cache_path: typing.Optional[pathlib.Path] = None

        self.cluster_metadata = inputs[0]

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
