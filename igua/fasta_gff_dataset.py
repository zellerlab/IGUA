import gc
import gzip
import io
import pathlib
import pickle
import re
import tarfile
import traceback
import typing
import uuid
from dataclasses import asdict, dataclass

import pandas as pd
import polars as pl
import rich.progress
from rich.console import Console

from .base_dataset import BaseDataset


_GZIP_MAGIC = b"\x1f\x8b"


class TarCache:
    """Cache for extracted tar members to avoid repeated extraction."""

    def __init__(self):
        self._cache: typing.Dict[str, bytes] = {}
        self._tar_handles: typing.Dict[pathlib.Path, tarfile.TarFile] = {}

    def get_member(self, tar_path: pathlib.Path, member_path: str) -> io.BytesIO:
        """Get tar member, using cache if available."""
        cache_key = f"{tar_path}::{member_path}"

        if cache_key not in self._cache:
            if tar_path not in self._tar_handles:
                self._tar_handles[tar_path] = tarfile.open(tar_path, "r:*")

            tf = self._tar_handles[tar_path]
            member = tf.extractfile(member_path)
            if member is None:
                raise FileNotFoundError(f"Member {member_path} not found in {tar_path}")

            self._cache[cache_key] = member.read()
            member.close()

        return io.BytesIO(self._cache[cache_key])

    def cleanup(self):
        """Close all tar handles and clear cache."""
        for tf in self._tar_handles.values():
            tf.close()
        self._tar_handles.clear()
        self._cache.clear()


def _parse_tar_path(
    path: pathlib.Path,
) -> typing.Union[tuple[pathlib.Path, str], tuple[None, None]]:
    """Parse a path that might point inside a tar archive.

    Args:
        path: Path that might contain a tar archive reference.

    Returns:
        Tuple of (tar_path, member_path) if tar archive found, else (None, None).
    """
    path_str = str(path)

    for ext in [".tar.gz", ".tar.bz2", ".tgz", ".tar"]:
        if ext in path_str:
            parts = path_str.split(ext, 1)
            tar_path = pathlib.Path(parts[0] + ext)
            if tar_path.exists():
                member_path = parts[1].lstrip("/")
                return tar_path, member_path

    return None, None


def smart_open(
    path: pathlib.Path, mode: str = "rb", tar_cache: typing.Optional[TarCache] = None
) -> typing.BinaryIO:
    """Open a file, handling both regular files and files inside tar archives.

    Uses caching for tar members to avoid repeated extraction.
    Supports gzip compression for both regular files and tar members.

    Args:
        path: Path to file (may be inside tar archive).
        mode: File mode ('rb' or 'rt').
        tar_cache: Optional TarCache instance for caching tar members.

    Returns:
        File-like object (binary mode).
    """
    tar_path, member_path = _parse_tar_path(path)

    if tar_path and member_path:
        if tar_cache is None:
            tar_cache = TarCache()
        member_data = tar_cache.get_member(tar_path, member_path)
        reader = io.BufferedReader(member_data)

        if reader.peek().startswith(_GZIP_MAGIC):
            reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore

        return reader  # type: ignore
    else:
        reader = open(path, "rb")
        if reader.peek().startswith(_GZIP_MAGIC):
            reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
        return reader  # type: ignore


@dataclass
class SystemCoordinates:
    """Genome coordinates for a gene cluster.

    Attributes:
        cluster_id: Cluster identifier.
        seq_id: Sequence/contig identifier.
        start_coord: Start coordinate on the sequence (1-based).
        end_coord: End coordinate on the sequence (1-based, inclusive).
        strand: Strand orientation ('+', '-', or '.').
        genes: List of gene identifiers in the system.
        fasta_file: Path to the genome FASTA file.
        valid: Whether the coordinates are valid.
        error_msg: Error message if coordinates are invalid.
    """

    cluster_id: str
    seq_id: str
    start_coord: int
    end_coord: int
    strand: str
    genes: typing.List[str]
    fasta_file: str
    valid: bool = True
    error_msg: typing.Optional[str] = None

    def to_dict(self) -> typing.Dict:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary representation of the coordinates.
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, data: typing.Dict) -> "SystemCoordinates":
        """Create from dictionary.

        Args:
            data: Dictionary containing coordinate data.

        Returns:
            SystemCoordinates instance.
        """
        return cls(**data)


class GFFIndex:
    """Fast, in-memory GFF index for gene feature lookup.

    Attributes:
        path: Path to the GFF file.
    """

    def __init__(
        self, gff_path: pathlib.Path, tar_cache: typing.Optional[TarCache] = None
    ):
        """Initialize GFF index.

        Args:
            gff_path: Path to the GFF file.
            tar_cache: Optional TarCache instance for caching tar members.
        """
        self.path = gff_path
        self._tar_cache = tar_cache
        self._index: typing.Dict[str, typing.Dict] = {}
        self._build_index()

    def _build_index(self):
        """Build comprehensive ID index from GFF file.

        Creates multiple lookup variants for each gene/CDS feature including
        ID, locus_tag, Name, gene, old_locus_tag, and protein_id attributes.
        Also creates prefixed variants (gene-, cds-) and underscore/tilde variants.
        """
        with smart_open(self.path, tar_cache=self._tar_cache) as reader:
            for line in io.TextIOWrapper(reader, encoding="utf-8"):
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] not in ["gene", "CDS"]:
                    continue

                seqid, _, ftype, start, end, _, strand, _, attrs = parts
                attr_dict = dict(
                    item.split("=", 1) for item in attrs.split(";") if "=" in item
                )

                feature = {
                    "seqid": seqid,
                    "type": ftype,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attributes": attr_dict,
                }

                for key in [
                    "ID",
                    "locus_tag",
                    "Name",
                    "gene",
                    "old_locus_tag",
                    "protein_id",
                ]:
                    if val := attr_dict.get(key):
                        for variant in [
                            val,
                            f"gene-{val}",
                            f"cds-{val}",
                            val.replace("_", "~"),
                            val.replace("~", "_"),
                        ]:
                            self._index[variant] = feature

    def get(self, gene_id: str) -> typing.Optional[typing.Dict]:
        """Get feature by gene ID with O(1) lookup.

        Args:
            gene_id: Gene identifier to look up.

        Returns:
            Feature dictionary if found, None otherwise.
        """
        return self._index.get(gene_id)

    def __contains__(self, gene_id: str) -> bool:
        """Check if gene ID exists in index.

        Args:
            gene_id: Gene identifier to check.

        Returns:
            True if gene ID exists, False otherwise.
        """
        return gene_id in self._index


def read_fasta(
    file_path: pathlib.Path, tar_cache: typing.Optional[TarCache] = None
) -> typing.Iterable[typing.Tuple[str, str, str]]:
    """Stream FASTA records from file.

    Handles both plain and gzip-compressed FASTA files.

    Args:
        file_path: Path to FASTA file (.fasta, .fa, .fna, or .gz).
        tar_cache: Optional TarCache instance for caching tar members.

    Yields:
        Tuple of (sequence_id, full_header, sequence_string).
    """
    with smart_open(file_path, tar_cache=tar_cache) as reader:
        name = None
        full_header = None
        sequence = []

        for line in io.TextIOWrapper(reader, encoding="utf-8"):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    yield name, full_header, "".join(sequence)
                full_header = line[1:]
                name = full_header.split()[0]
                sequence = []
            else:
                sequence.append(line)

        if name is not None:
            yield name, full_header, "".join(sequence)


class ProteinIndex:
    """Lazy-loading protein sequence index for efficient lookup.

    Only loads proteins as needed, suitable for large FASTA files.

    Attributes:
        path: Path to the protein FASTA file.
    """

    def __init__(
        self, protein_fasta: pathlib.Path, tar_cache: typing.Optional[TarCache] = None
    ):
        """Initialize protein index without loading sequences.

        Args:
            protein_fasta: Path to protein FASTA file.
            tar_cache: Optional TarCache instance for caching tar members.
        """
        self.path = protein_fasta
        self._tar_cache = tar_cache
        self._sequences: typing.Dict[str, str] = {}
        self._headers: typing.Dict[str, str] = {}
        self._loaded = False

    def _ensure_loaded(self, gene_ids: typing.Optional[set] = None):
        """Load protein sequences on demand.

        Args:
            gene_ids: Optional set of gene IDs to load. If None, loads all.
        """
        if self._loaded:
            return

        for seq_id, full_header, sequence in read_fasta(
            self.path, tar_cache=self._tar_cache
        ):
            if gene_ids is None or seq_id in gene_ids:
                self._sequences[seq_id] = sequence
                self._headers[seq_id] = full_header

        self._loaded = True

    def load_subset(self, gene_ids: set):
        """Pre-load only specific proteins by ID.

        Args:
            gene_ids: Set of gene IDs to load.
        """
        self._ensure_loaded(gene_ids)

    def get(self, protein_id: str) -> typing.Optional[str]:
        """Get protein sequence by ID.

        Args:
            protein_id: Protein identifier.

        Returns:
            Protein sequence if found, None otherwise.
        """
        if not self._loaded:
            self._ensure_loaded()
        return self._sequences.get(protein_id)

    def get_with_fallback(self, gene_id: str) -> typing.Optional[str]:
        """Get protein sequence with fallback header search."""
        seq = self.get(gene_id)
        if seq:
            return seq

        opener = gzip.open if self.path.suffix == ".gz" else open

        with opener(self.path, "rt") as f:
            seq_id = None
            sequence = []
            full_header = None

            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if seq_id and full_header:
                        for attr in ["locus_tag", "ID", "Name", "gene"]:
                            if re.search(
                                rf"\[{attr}=({re.escape(gene_id)})\]", full_header
                            ):
                                return "".join(sequence)

                    full_header = line[1:]
                    seq_id = full_header.split()[0]
                    sequence = []
                else:
                    sequence.append(line)

            if seq_id and full_header:
                for attr in ["locus_tag", "ID", "Name", "gene"]:
                    if re.search(rf"\[{attr}=({re.escape(gene_id)})\]", full_header):
                        return "".join(sequence)

        return None


class GenomeContext:
    """Immutable data container for genome/MAG file paths and metadata."""

    def __init__(
        self,
        genome_id: typing.Optional[str],
        cluster_tsv: pathlib.Path,
        gff_file: pathlib.Path,
        genome_fasta: pathlib.Path,
        protein_fasta: pathlib.Path,
        column_mapping: typing.Dict[str, str],
        system_loader: typing.Callable[
            [pathlib.Path, Console], pl.DataFrame
        ],  # Add this
    ):

        self.genome_id = genome_id if genome_id else str(uuid.uuid4())[:8]
        self.cluster_tsv = pathlib.Path(cluster_tsv)
        self.gff_file = pathlib.Path(gff_file)
        self.genome_fasta = pathlib.Path(genome_fasta)
        self.protein_fasta = pathlib.Path(protein_fasta)
        self.column_mapping = column_mapping
        self.system_loader = system_loader

        self.missing_files = [
            f"{name}: {path}"
            for path, name in [
                (self.cluster_tsv, "cluster_tsv"),
                (self.gff_file, "gff_file"),
                (self.genome_fasta, "genome_fasta"),
                (self.protein_fasta, "protein_fasta"),
            ]
            if not path.exists() and "tar.gz" not in str(path)
        ]

    def __repr__(self):
        return (
            f"<GenomeContext "
            f"genome_id={self.genome_id!r} "
            f"files={4 - len(self.missing_files)}/4>"
        )

    def is_valid(self) -> bool:
        return len(self.missing_files) == 0


class GenomeResources:
    """Manages lazy-loading and caching of genome resources.

    Format-agnostic - uses adapter for format-specific operations.

    Attributes:
        context: Genome context with file paths and adapter.
        console: Rich console for logging.
    """

    def __init__(self, context: GenomeContext, console: Console):
        """Initialize genome resources manager.

        Args:
            context: Genome context with file paths and adapter.
            console: Rich console for logging.
        """
        self.context = context
        self.console = console
        self._cluster_df: typing.Optional[pl.DataFrame] = None
        self._protein_idx: typing.Optional[ProteinIndex] = None
        self._gff_db: typing.Optional[GFFIndex] = None
        self._coordinates_cache: typing.Optional[typing.List[SystemCoordinates]] = None
        self._tar_cache = TarCache()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._tar_cache.cleanup()
        return False

    @property
    def cluster_df(self) -> pl.DataFrame:
        """Load and filter clusters TSV (lazy-loaded, cached).

        Delegates format-specific logic to adapter.

        Returns:
            Polars DataFrame with filtered systems.
        """
        if self._cluster_df is not None:
            return self._cluster_df

        df = self.context.system_loader(self.context.cluster_tsv, self.console)

        col_map = self.context.column_mapping
        cluster_col = col_map["cluster_id"]

        if (
            n_dup := df.filter(pl.col(cluster_col).is_duplicated())
            .select(cluster_col)
            .unique()
            .height
        ) > 0:
            dup_ids = (
                df.filter(pl.col(cluster_col).is_duplicated())
                .select(cluster_col)
                .unique()
                .to_series()
                .to_list()
            )
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] {n_dup} duplicate cluster/s in [bold cyan]{self.context.genome_id}[/]: "
                f"[cyan]{', '.join(dup_ids[:5])}{'...' if n_dup > 5 else ''}[/]"
            )
            df = df.unique(subset=[cluster_col], keep="first")

        self._cluster_df = df
        return self._cluster_df

    @property
    def protein_idx(self) -> ProteinIndex:
        """Get protein index (lazy-loaded, cached).

        Returns:
            ProteinIndex instance.
        """
        if self._protein_idx is None:
            self._protein_idx = ProteinIndex(
                self.context.protein_fasta, tar_cache=self._tar_cache
            )
        return self._protein_idx

    @property
    def gff_db(self) -> GFFIndex:
        """Get GFF index (lazy-loaded, cached).

        Returns:
            GFFIndex instance.
        """
        if self._gff_db is None:
            self._gff_db = GFFIndex(self.context.gff_file, tar_cache=self._tar_cache)
        return self._gff_db

    @property
    def coordinates(self) -> typing.List[SystemCoordinates]:
        """Build and cache system coordinates (lazy-loaded, cached).

        Returns:
            List of SystemCoordinates for all systems.
        """
        if self._coordinates_cache is None:
            self._coordinates_cache = self._build_coordinates()
        return self._coordinates_cache

    def _build_coordinates(self) -> typing.List[SystemCoordinates]:
        """Build coordinates for all clusters.

        Returns:
            List of SystemCoordinates.
        """
        coordinates = []
        for row in self.cluster_df.iter_rows(named=True):
            coord = self._parse_system_coordinates(row)
            coordinates.append(coord)
        return coordinates

    def _parse_system_coordinates(self, row: dict) -> SystemCoordinates:
        """Parse coordinates for a single system.

        Uses adapter's column mapping for format flexibility.

        Args:
            row: Dictionary containing system data from TSV row.

        Returns:
            SystemCoordinates instance.
        """
        col_map = self.context.column_mapping
        cluster_id = row[col_map["cluster_id"]]

        sys_beg_gene = row.get(col_map["start_gene"])
        sys_end_gene = row.get(col_map["end_gene"])

        if sys_beg_gene is None or sys_end_gene is None:
            return self._invalid_coord(
                cluster_id,
                [],
                f"Missing '{col_map['sys_beg']}' or '{col_map['sys_end']}' columns in systems TSV",
            )

        try:
            gene_list = [
                g.strip()
                for g in str(row[col_map["genes_in_cluster"]]).split(",")
                if g.strip()
            ]
        except (KeyError, AttributeError):
            return self._invalid_coord(
                cluster_id,
                [],
                f"Missing '{col_map['genes_in_cluster']}' column",
            )

        if not gene_list:
            return self._invalid_coord(cluster_id, [], "Empty gene list")

        beg_feature = self.gff_db.get(sys_beg_gene)
        end_feature = self.gff_db.get(sys_end_gene)

        if beg_feature is None:
            return self._invalid_coord(
                cluster_id, gene_list, f"Start gene '{sys_beg_gene}' not found in GFF"
            )

        if end_feature is None:
            return self._invalid_coord(
                cluster_id, gene_list, f"End gene '{sys_end_gene}' not found in GFF"
            )

        seq_id_beg = beg_feature["seqid"]
        seq_id_end = end_feature["seqid"]

        if seq_id_beg != seq_id_end:
            return self._invalid_coord(
                cluster_id,
                gene_list,
                f"Start and end genes on different contigs/seqs: {seq_id_beg} vs {seq_id_end} for [bold cyan]{str(self.context.genome_id)}[/]",
            )

        seq_id = seq_id_beg

        start = min(
            beg_feature["start"],
            beg_feature["end"],
            end_feature["start"],
            end_feature["end"],
        )
        end = max(
            beg_feature["start"],
            beg_feature["end"],
            end_feature["start"],
            end_feature["end"],
        )

        strand = beg_feature["strand"]

        region_size = end - start + 1
        if region_size > 1e5:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] Cluster [bold cyan]{cluster_id}[/] unusually large: {region_size} bp"
            )
        elif region_size < 50:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] Cluster [bold cyan]{cluster_id}[/] unusually small: {region_size} bp"
            )

        return SystemCoordinates(
            cluster_id=cluster_id,
            seq_id=seq_id,
            start_coord=start,
            end_coord=end,
            strand=strand,
            genes=gene_list,
            fasta_file=str(self.context.genome_fasta),
            valid=True,
        )

    def _invalid_coord(
        self, cluster_id: str, genes: typing.List[str], error: str, seq_id: str = ""
    ) -> SystemCoordinates:
        """Create an invalid SystemCoordinates object.

        Args:
            cluster_id: Cluster identifier.
            genes: List of gene identifiers.
            error: Error message describing the issue.
            seq_id: Sequence identifier (default: empty string).

        Returns:
            Invalid SystemCoordinates instance with error message.
        """
        self.console.print(
            f"[bold yellow]{'Warning':>12}[/] {error} for cluster [bold cyan]{cluster_id}[/]"
        )
        return SystemCoordinates(
            cluster_id=cluster_id,
            seq_id=seq_id,
            start_coord=0,
            end_coord=0,
            strand="",
            genes=genes,
            fasta_file=str(self.context.genome_fasta),
            valid=False,
            error_msg=error,
        )

    def extract_genome_sequences(
        self,
        output_file: typing.TextIO,
        verbose: bool = False,
        progress_task: typing.Optional[int] = None,
        progress: typing.Optional[rich.progress.Progress] = None,
    ) -> typing.List[typing.Tuple[str, int, str]]:
        """Extract genome sequences by streaming FASTA file."""
        coordinates = self.coordinates
        valid_coords = [c for c in coordinates if c.valid]

        if not valid_coords:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] No valid clusters to extract for [bold cyan]{self.context.genome_id}[/]"
            )
            return []

        contig_groups: typing.Dict[str, typing.List[SystemCoordinates]] = {}
        for coord in valid_coords:
            contig_groups.setdefault(coord.seq_id, []).append(coord)

        num_contigs = len(contig_groups)

        if progress and progress_task is not None:
            progress.update(
                progress_task,
                description=f"[bold blue]{'Processing':>9}[/] {len(valid_coords)} clusters across {num_contigs} contigs/seqs for [bold cyan]{str(self.context.genome_id)}[/]",
            )
        else:
            self.console.print(
                f"[bold blue]{'Processing':>12}[/] {len(valid_coords)} clusters across {num_contigs} contigs/seqs for [bold cyan]{str(self.context.genome_id)}[/]"
            )

        results = []

        for seq_id, _, sequence in read_fasta(
            self.context.genome_fasta, tar_cache=self._tar_cache
        ):
            if seq_id not in contig_groups:
                continue

            if verbose:
                self.console.print(
                    f"[bold blue]{'Loading':>12}[/] contig {seq_id} ({len(contig_groups[seq_id])} clusters) for [bold cyan]{self.context.genome_id}[/]"
                )

            for coord in contig_groups[seq_id]:
                try:
                    subseq = sequence[coord.start_coord - 1 : coord.end_coord]
                    output_file.write(f">{coord.cluster_id}\n{subseq}\n")

                    if verbose:
                        self.console.print(
                            f"[bold blue]{'Extracted':>12}[/] [cyan]{coord.cluster_id}[/] ({len(subseq)} bp) for [bold cyan]{self.context.genome_id}[/]"
                        )

                    results.append((coord.cluster_id, len(subseq), coord.fasta_file))
                except Exception as e:
                    self.console.print(
                        f"[bold red]{'Error':>12}[/] Failed to extract [cyan]{coord.cluster_id}[/] for [bold cyan]{self.context.genome_id}[/]: {e}"
                    )

            del contig_groups[seq_id]

            if not contig_groups:
                break

        if contig_groups:
            for seq_id in contig_groups:
                self.console.print(
                    f"[bold red]{'Error':>12}[/] Contig {seq_id} not found in genome [bold cyan]{self.context.genome_id}[/]"
                )

        return results

    def extract_proteins_from_coordinates(
        self,
        coordinates: typing.List[SystemCoordinates],
        output_file: typing.TextIO,
        verbose: bool = False,
        progress_task: typing.Optional[int] = None,
        progress: typing.Optional[rich.progress.Progress] = None,
    ) -> typing.Dict[str, int]:
        """Extract protein sequences from coordinates."""
        valid_coords = [c for c in coordinates if c.valid]

        if not valid_coords:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] No valid clusters for [bold cyan]{self.context.genome_id}[/]"
            )
            return {}

        all_gene_ids = set()
        for coord in valid_coords:
            all_gene_ids.update(coord.genes)

        self.protein_idx.load_subset(all_gene_ids)

        total_genes = len(all_gene_ids)

        if progress and progress_task is not None:
            progress.update(
                progress_task,
                description=f"[bold blue]{'Processing':>9}[/] {total_genes} proteins from {len(valid_coords)} clusters for [bold cyan]{str(self.context.genome_id)}[/]",
            )
        else:
            self.console.print(
                f"[bold blue]{'Processing':>12}[/] {total_genes} proteins from {len(valid_coords)} clusters for [bold cyan]{str(self.context.genome_id)}[/]"
            )

        protein_sizes = {}
        for coord in valid_coords:
            for gene_id in coord.genes:
                if seq := self.protein_idx.get_with_fallback(gene_id):
                    protein_id = f"{coord.cluster_id}__{gene_id}"
                    output_file.write(f">{protein_id}\n{seq}\n")
                    protein_sizes[protein_id] = len(seq)
                else:
                    if verbose:
                        self.console.print(
                            f"[bold yellow]{'Warning':>12}[/] Protein {gene_id} not found for cluster [bold cyan]{coord.cluster_id}[/]"
                        )

        self.console.print(
            f"[bold blue]{'Extracted':>12}[/] {len(protein_sizes)} proteins from "
            f"{len(valid_coords)} clusters for [bold cyan]{self.context.genome_id}[/]"
        )

        return protein_sizes


class ClusterMetadataCache:
    """Manages caching of cluster metadata to avoid redundant processing.

    Uses batched pickle format for memory-efficient serialization/deserialization.

    Attributes:
        cache_path: Path to cache file.
        batch_size: Number of genomes per batch (default 500).
    """

    def __init__(self, cache_path: pathlib.Path, batch_size: int = 500):
        """Initialize metadata cache.

        Args:
            cache_path: Path to cache file.
            batch_size: Number of genomes to batch together.
            _metadata: Metadata cache.
        """
        self.cache_path = cache_path
        self.batch_size = batch_size
        self._metadata: typing.Optional[typing.Dict] = None

    def save_streaming(self, genomes_iterator, total: int) -> None:
        """Save genomes in batches to reduce memory usage.

        Args:
            genomes_iterator: Iterator yielding genome metadata dictionaries.
            total: Total number of genomes (for compression decision).
        """
        batches = []
        batch = []

        for genome in genomes_iterator:
            batch.append(genome)
            if len(batch) >= self.batch_size:
                batches.append(batch)
                batch = []

        if batch:
            batches.append(batch)

        use_compression = total > 10000
        cache_path = (
            self.cache_path.with_suffix(".pkl.gz")
            if use_compression
            else self.cache_path
        )

        open_fn = gzip.open if use_compression else open
        with open_fn(cache_path, "wb") as f:
            pickle.dump(
                {
                    "genomes": batches,
                    "batch_size": self.batch_size,
                    "compressed": use_compression,
                },
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

        if use_compression:
            self.cache_path = cache_path

    def iter_genomes(self):
        """Stream genomes batch-by-batch without loading all into memory.

        Yields:
            Individual genome metadata dictionaries.
        """
        compressed_path = self.cache_path.with_suffix(".pkl.gz")
        use_compression = compressed_path.exists()
        cache_path = compressed_path if use_compression else self.cache_path

        open_fn = gzip.open if use_compression else open
        with open_fn(cache_path, "rb") as f:
            data = pickle.load(f)

            if isinstance(data.get("genomes"), list) and data.get("batch_size"):
                for batch in data["genomes"]:
                    for genome in batch:
                        yield genome
                    del batch
                    gc.collect()
            else:
                for genome in data.get("genomes", []):
                    yield genome

    def exists(self) -> bool:
        """Check if cache file exists.

        Returns:
            True if cache file exists, False otherwise.
        """
        return (
            self.cache_path.exists() or self.cache_path.with_suffix(".pkl.gz").exists()
        )

    def clear(self) -> None:
        """Clear cache file and memory."""
        if self.cache_path.exists():
            self.cache_path.unlink()
        compressed_path = self.cache_path.with_suffix(".pkl.gz")
        if compressed_path.exists():
            compressed_path.unlink()
        self._metadata = None


class FastaGFFDataset(BaseDataset):
    """FastaGFF dataset class."""

    def __init__(
        self,
        inputs: typing.List[pathlib.Path],
        column_mapping: typing.Optional[typing.Dict[str, str]] = None,
    ) -> None:
        """Initialize the FastaGFFDataset class.

        Args:
            inputs: List of input paths (metadata TSV or individual files).
            column_mapping: Custom column mapping. If None, uses default generic mapping.
        """
        super().__init__(inputs)

        self.column_mapping = column_mapping or {
            "cluster_id": "cluster_id",
            "start_gene": "start_gene",
            "end_gene": "end_gene",
            "genes_in_cluster": "genes_in_cluster",
        }

        self.verbose: bool = False
        self.gff_cache_dir: typing.Optional[pathlib.Path] = None
        self._metadata_cache_path: typing.Optional[pathlib.Path] = None
        self.cluster_metadata = inputs[0]

    def _create_genome_context(self, row: typing.Dict, genome_id: str) -> GenomeContext:
        return GenomeContext(
            genome_id=genome_id,
            cluster_tsv=pathlib.Path(row["cluster_tsv"]),
            gff_file=pathlib.Path(row["gff_file"]),
            genome_fasta=pathlib.Path(row["genome_fasta_file"]),
            protein_fasta=pathlib.Path(row["protein_fasta_file"]),
            column_mapping=self.column_mapping,
            system_loader=self._load_and_filter_systems,
        )

    def _load_and_filter_systems(
        self, tsv_path: pathlib.Path, console: Console
    ) -> pl.DataFrame:
        """Load systems TSV without filtering (generic behavior)."""
        df = pl.read_csv(tsv_path, separator="\t")
        return df

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
    ) -> pd.DataFrame:
        """Extracts nucleotide sequences from gene clusters and caches metadata."""
        progress.console.print(
            f"[bold blue]{'Using':>12}[/] cluster metadata file: [magenta]{self.cluster_metadata}[/]"
        )

        df = pl.read_csv(self.cluster_metadata, separator="\t")
        self._metadata_cache_path = output.parent / f".{output.stem}_metadata.json"

        results = []

        def metadata_generator():
            """Generator that yields metadata for each genome."""
            genome_count = 0
            task = progress.add_task(
                f"[bold blue]{'Processing':>9}[/] gene clusters", total=len(df)
            )

            build_task = None
            validate_task = None
            processing_task = None

            with open(output, "w") as dst:
                for row in df.iter_rows(named=True):
                    genome_id = row.get("genome_id") or f"genome_{genome_count:07}"
                    genome_count += 1
                    progress.update(
                        task,
                        description=f"[bold blue]{'Processing':>9}[/] genome: [bold cyan]{genome_id}",
                    )

                    context = self._create_genome_context(row, genome_id)

                    if not context.is_valid():
                        progress.console.print(
                            f"[bold yellow]{'Missing':>12}[/] files for [bold cyan]{genome_id}[/]"
                        )
                        progress.update(task, advance=1)
                        continue

                    try:
                        with GenomeResources(context, progress.console) as resources:
                            build_task = progress.add_task(
                                f"[bold blue]{'Building':>9}[/] cluster coordinates for [bold cyan]{str(genome_id)}[/]",
                                total=None,
                            )

                            coordinates = resources.coordinates
                            progress.remove_task(build_task)
                            build_task = None

                            valid_count = sum(1 for c in coordinates if c.valid)
                            
                            validate_task = progress.add_task(
                                f"[bold green]{'Validated':>9}[/] {valid_count}/{len(coordinates)} clusters for [bold cyan]{str(genome_id)}[/]",
                                total=None,
                            )

                            progress.remove_task(validate_task)
                            validate_task = None

                            processing_task = progress.add_task(
                                f"[bold blue]{'Extracting':>9}[/] cluster genome sequences for [bold cyan]{str(genome_id)}[/]",
                                total=None,
                            )

                            extraction_results = resources.extract_genome_sequences(
                                dst,
                                verbose=self.verbose,
                                progress_task=processing_task,
                                progress=progress,
                            )

                            progress.remove_task(processing_task)
                            processing_task = None

                            contig_groups = {}
                            for coord in coordinates:
                                if coord.valid:
                                    contig_groups.setdefault(coord.seq_id, []).append(
                                        coord
                                    )
                            num_contigs = len(contig_groups)

                            progress.console.print(
                                f"[bold blue]{'Extracted':>12}[/] {len(extraction_results)} gene clusters across {num_contigs} contigs/seqs for [bold cyan]{context.genome_id}[/]"
                            )

                            for coord in coordinates:
                                if coord.valid:
                                    results.append(
                                        (
                                            coord.cluster_id,
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
                        if build_task is not None:
                            progress.remove_task(build_task)
                            build_task = None
                        if validate_task is not None:
                            progress.remove_task(validate_task)
                            validate_task = None
                        if processing_task is not None:
                            progress.remove_task(processing_task)
                            processing_task = None

                        progress.console.print(
                            f"[bold red]{'Error':>12}[/] processing [bold cyan]{genome_id}[/]: {e}"
                        )
                        traceback.print_exc()

                    progress.update(task, advance=1)

            progress.remove_task(task)

        cache = ClusterMetadataCache(self._metadata_cache_path)
        cache.save_streaming(metadata_generator(), len(df))

        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(results):,} gene clusters from {len(df):,} genomes/MAGs"
        )
        progress.console.print(
            f"[bold blue]{'Cached':>12}[/] metadata to [magenta]{self._metadata_cache_path.name}[/]"
        )

        return pd.DataFrame(
            data=results, columns=["cluster_id", "cluster_length", "filename"]
        ).set_index("cluster_id")

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Extracts protein sequences using cached metadata."""
        if self._metadata_cache_path and self._metadata_cache_path.exists():
            progress.console.print(
                f"[bold blue]{'Using':>12}[/] cached metadata from sequence extraction"
            )
            return self._extract_proteins_from_cache(progress, output, representatives)

        progress.console.print(
            f"[bold yellow]{'Warning':>12}[/] No cached metadata found, reading cluster metadata again"
        )
        progress.console.print(
            f"[bold blue]{'Using':>12}[/] cluster metadata file: [magenta]{self.cluster_metadata}[/]"
        )

        try:
            df = pl.read_csv(self.cluster_metadata, separator="\t")
            return self._extract_proteins_direct(progress, df, output, representatives)
        except Exception as e:
            progress.console.print(
                f"[bold red]{'Error':>12}[/] reading cluster metadata: {e}"
            )
            return {}

    def _extract_proteins_from_cache(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
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
                    description=f"[bold blue]{'Processing':>9}[/] genome: [bold cyan]{genome_id}",
                )

                try:
                    proteins = self._extract_proteins_from_metadata(
                        genome_metadata,
                        dst,
                        progress.console,
                        representatives,
                        progress,
                    )
                    protein_sizes.update(proteins)
                except Exception as e:
                    progress.console.print(
                        f"[bold red]{'Error':>12}[/] processing [bold cyan]{genome_id}[/]: {e}"
                    )
                    traceback.print_exc()

                progress.update(task, advance=1)

            progress.remove_task(task)

        self._log_protein_summary(progress, protein_sizes, representatives)
        return protein_sizes

    def _extract_proteins_direct(
        self,
        progress: rich.progress.Progress,
        df: pl.DataFrame,
        output: pathlib.Path,
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
                    description=f"[bold blue]{'Processing':>9}[/] genome: [bold cyan]{genome_id}",
                )

                context = self._create_genome_context(row, genome_id)

                if not context.is_valid():
                    progress.console.print(
                        f"[bold yellow]{'Missing':>12}[/] files for {genome_id}"
                    )
                    progress.update(task, advance=1)
                    continue

                try:
                    with GenomeResources(context, progress.console) as resources:
                        coordinates = resources.coordinates

                        if representatives:
                            rep_set = (
                                set(representatives)
                                if not isinstance(representatives, set)
                                else representatives
                            )
                            coordinates = [
                                c for c in coordinates if c.cluster_id in rep_set
                            ]

                        proteins = resources.extract_proteins_from_coordinates(
                            coordinates, dst, verbose=self.verbose
                        )
                        protein_sizes.update(proteins)

                except Exception as e:
                    progress.console.print(
                        f"[bold red]{'Error':>12}[/] processing {genome_id}: {e}"
                    )
                    traceback.print_exc()

                progress.update(task, advance=1)

            progress.remove_task(task)

        self._log_protein_summary(progress, protein_sizes, representatives)
        return protein_sizes

    def _extract_proteins_from_metadata(
        self,
        metadata: typing.Dict,
        output_file: typing.TextIO,
        console: Console,
        representatives: typing.Optional[typing.Container[str]] = None,
        progress: typing.Optional[rich.progress.Progress] = None,
    ) -> typing.Dict[str, int]:
        """Extract protein sequences using pre-computed metadata."""
        genome_id = metadata["genome_id"]
        protein_fasta = pathlib.Path(metadata["protein_fasta"])
        coordinates = [SystemCoordinates.from_dict(c) for c in metadata["coordinates"]]

        if representatives:
            rep_set = (
                set(representatives)
                if not isinstance(representatives, set)
                else representatives
            )
            coordinates = [c for c in coordinates if c.cluster_id in rep_set]

        valid_coords = [c for c in coordinates if c.valid]

        if not valid_coords:
            console.print(
                f"[bold yellow]{'Warning':>12}[/] No clusters to extract for [bold cyan]{str(genome_id)}[/]"
            )
            return {}

        tar_cache = TarCache()
        processing_task = None

        try:
            all_gene_ids = set()
            for coord in valid_coords:
                all_gene_ids.update(coord.genes)

            protein_idx = ProteinIndex(protein_fasta, tar_cache=tar_cache)
            protein_idx.load_subset(all_gene_ids)

            total_genes = len(all_gene_ids)

            if progress:
                processing_task = progress.add_task(
                    f"[bold blue]{'Processing':>9}[/] {total_genes} proteins from {len(valid_coords)} clusters for [bold cyan]{str(genome_id)}[/]",
                    total=None,
                )
            else:
                console.print(
                    f"[bold blue]{'Processing':>12}[/] {total_genes} proteins from {len(valid_coords)} clusters for [bold cyan]{str(genome_id)}[/]"
                )

            protein_sizes = {}
            for coord in valid_coords:
                for gene_id in coord.genes:
                    if seq := protein_idx.get_with_fallback(gene_id):
                        protein_id = f"{coord.cluster_id}__{gene_id}"
                        output_file.write(f">{protein_id}\n{seq}\n")
                        protein_sizes[protein_id] = len(seq)
                    else:
                        if self.verbose:
                            console.print(
                                f"[bold yellow]{'Warning':>12}[/] Protein {gene_id} not found for cluster [bold cyan]{coord.cluster_id}[/]"
                            )

            if processing_task is not None:
                progress.remove_task(processing_task)
                processing_task = None

            console.print(
                f"[bold blue]{'Extracted':>12}[/] {len(protein_sizes)} proteins from "
                f"{len(valid_coords)} clusters for [bold cyan]{str(genome_id)}[/]"
            )

            return protein_sizes

        except Exception as e:
            if processing_task is not None:
                progress.remove_task(processing_task)
                processing_task = None

            console.print(
                f"[bold red]{'Error':>12}[/] Fatal error for {genome_id}: {e}"
            )
            traceback.print_exc()
            return {}
        finally:
            tar_cache.cleanup()

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
