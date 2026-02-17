import pathlib
import typing

import polars as pl
from rich.console import Console

from .fasta_gff import FastaGFFDataset


class DefenseFinderDataset(FastaGFFDataset):
    """DefenseFinder-specific dataset class."""

    def __init__(
        self,
        inputs: typing.List[pathlib.Path],
        activity_filter: str = "all",
    ) -> None:
        """Initialize DefenseFinder dataset.

        Args:
            inputs: List of input file paths (typically metadata TSV).
            activity_filter: Filter by activity type. Options: 'all',
                'defense', 'antidefense'.
        """
        defensefinder_mapping = {
            "cluster_id": "sys_id",
            "start_gene": "sys_beg",
            "end_gene": "sys_end",
            "genes_in_cluster": "protein_in_syst",
        }

        super().__init__(inputs, column_mapping=defensefinder_mapping)
        self.activity_filter = activity_filter

    def _load_and_filter_systems(
        self, tsv_path: pathlib.Path, console: Console
    ) -> pl.DataFrame:
        """Load DefenseFinder systems with activity filtering.
        
        Args:
            tsv_path: Path to the DefenseFinder systems TSV file.
            console: Rich console for logging output.
            
        Returns:
            Filtered Polars DataFrame containing system data.
        """
        df = pl.read_csv(tsv_path, separator="\t")
        original_count = len(df)

        if self.activity_filter.lower() != "all":
            if "activity" in df.columns:
                df = df.filter(
                    pl.col("activity").str.to_lowercase()
                    == self.activity_filter.lower()
                )
                console.print(
                    f"[bold green]{'Filtered':>12}[/] {original_count} → {len(df)} "
                    f"([bold cyan]{self.activity_filter}[/] only)"
                )
            else:
                console.print(
                    f"[bold yellow]{'Warning':>12}[/] No 'activity' column found"
                )

        return df