import collections
import gc
import gzip
import io
import pathlib
import pickle
import re
import tarfile
import typing
import uuid
from dataclasses import asdict, dataclass

import Bio.SeqIO
import pandas
import rich.progress
from rich.console import Console

from .base import BaseDataset, Cluster, Protein
from .._utils import zopen


class TarCache:
    """Cache for extracted tar members to avoid repeated extraction."""

    def __init__(self):
        """Initialize empty cache."""
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
    path: pathlib.Path,
    mode: str = "rb",
    tar_cache: typing.Optional[TarCache] = None,
) -> typing.BinaryIO:
    """Open file, handling regular files and tar archives.

    Uses caching for tar members to avoid repeated extraction.
    Supports gzip compression for both file types.

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
        return zopen(member_data)
    else:
        return zopen(path)


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

        Creates multiple lookup variants for each gene/CDS feature
        including ID, locus_tag, Name, gene, old_locus_tag, and
        protein_id attributes. Also creates prefixed variants
        (gene-, cds-) and underscore/tilde variants.
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
        with io.TextIOWrapper(reader, encoding="utf-8") as text_reader:
            yield from Bio.SeqIO.parse(text_reader, "fasta")


class ProteinIndex:
    """Lazy-loading protein sequence index for efficient lookup."""

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

        for record in read_fasta(self.path, tar_cache=self._tar_cache):
            if gene_ids is None or record.id in gene_ids:
                self._sequences[record.id] = str(record.seq)
                self._headers[record.id] = record.description

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
        """Get protein sequence with fallback to header search.

        Args:
            gene_id: Gene identifier to search for.

        Returns:
            Protein sequence if found, None otherwise.
        """
        seq = self.get(gene_id)
        if seq:
            return seq

        with smart_open(self.path) as reader:
            with io.TextIOWrapper(reader, encoding="utf-8") as text_reader:
                for record in Bio.SeqIO.parse(text_reader, "fasta"):
                    for attr in ["locus_tag", "ID", "Name", "gene"]:
                        if re.search(rf"\[{attr}=({re.escape(gene_id)})\]", record.description):
                            return str(record.seq)

        return None


class GenomeContext:
    """Data container for genome/MAG file paths and metadata."""

    def __init__(
        self,
        genome_id: typing.Optional[str],
        cluster_tsv: pathlib.Path,
        gff_file: pathlib.Path,
        genome_fasta: pathlib.Path,
        protein_fasta: pathlib.Path,
        column_mapping: typing.Dict[str, str],
        system_loader: typing.Callable[[pathlib.Path, Console], pandas.DataFrame],
    ):
        """Initialize genome context.

        Args:
            genome_id: Genome identifier (generates UUID if None).
            cluster_tsv: Path to cluster metadata TSV file.
            gff_file: Path to GFF annotation file.
            genome_fasta: Path to genome FASTA file.
            protein_fasta: Path to protein FASTA file.
            column_mapping: Column name mapping for TSV parsing.
            system_loader: Function to load and filter systems.
        """
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
        """Return string representation."""
        return (
            f"<GenomeContext "
            f"genome_id={self.genome_id!r} "
            f"files={4 - len(self.missing_files)}/4>"
        )

    def is_valid(self) -> bool:
        """Check if all required files exist.

        Returns:
            True if all files are present, False otherwise.
        """
        return len(self.missing_files) == 0


class GenomeResources:
    """Manages lazy-loading and caching of genome resources."""

    def __init__(self, context: GenomeContext, console: Console):
        """Initialize genome resources manager.

        Args:
            context: Genome context with file paths and adapter.
            console: Rich console for logging.
        """
        self.context = context
        self.console = console
        self._cluster_df: typing.Optional[pandas.DataFrame] = None
        self._protein_idx: typing.Optional[ProteinIndex] = None
        self._gff_db: typing.Optional[GFFIndex] = None
        self._coordinates_cache: typing.Optional[typing.List[SystemCoordinates]] = None
        self._tar_cache = TarCache()

    def __enter__(self):
        """Enter context manager."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager and cleanup resources."""
        self._tar_cache.cleanup()
        return False

    @property
    def cluster_df(self) -> pandas.DataFrame:
        """Load and filter clusters TSV.

        Returns:
            Polars DataFrame with filtered systems.
        """
        if self._cluster_df is not None:
            return self._cluster_df

        df = self.context.system_loader(self.context.cluster_tsv, self.console)

        col_map = self.context.column_mapping
        cluster_col = col_map["cluster_id"]

        dup = df.duplicated(cluster_col)
        n_dup = dup.sum()
        if n_dup > 0:
            dup_ids = df[dup][cluster_col]
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] {n_dup} duplicate cluster/s in [bold cyan]{self.context.genome_id}[/]: "
                f"[cyan]{', '.join(dup_ids[:5])}{'...' if n_dup > 5 else ''}[/]"
            )
        # if (
        #     n_dup := df.filter(pl.col(cluster_col).is_duplicated())
        #     .select(cluster_col)
        #     .unique()
        #     .height
        # ) > 0:
        #     dup_ids = (
        #         df.filter(pl.col(cluster_col).is_duplicated())
        #         .select(cluster_col)
        #         .unique()
        #         .to_series()
        #         .to_list()
        #     )
        #     self.console.print(
        #         f"[bold yellow]{'Warning':>12}[/] {n_dup} duplicate cluster/s in [bold cyan]{self.context.genome_id}[/]: "
        #         f"[cyan]{', '.join(dup_ids[:5])}{'...' if n_dup > 5 else ''}[/]"
        #     )
        #     df = df.unique(subset=[cluster_col], keep="first")
            df = df[~dup]

        self._cluster_df = df
        return self._cluster_df

    @property
    def protein_idx(self) -> ProteinIndex:
        """Get protein index.

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
        """Get GFF index.

        Returns:
            GFFIndex instance.
        """
        if self._gff_db is None:
            self._gff_db = GFFIndex(self.context.gff_file, tar_cache=self._tar_cache)
        return self._gff_db

    @property
    def coordinates(self) -> typing.List[SystemCoordinates]:
        """Build and cache system coordinates.

        Returns:
            List of SystemCoordinates for all systems.
        """
        if self._coordinates_cache is None:
            self._coordinates_cache = self._build_coordinates()
        return self._coordinates_cache

    def _build_coordinates(self) -> typing.List[SystemCoordinates]:
        """Parse coordinates from cluster DataFrame.

        Returns:
            List of SystemCoordinates for each cluster.
        """
        coordinates = []
        for row in self.cluster_df.itertuples():
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
        cluster_id = getattr(row, col_map['cluster_id'])

        sys_beg_gene = getattr(row, col_map['start_gene'])
        sys_end_gene = getattr(row, col_map['end_gene'])

        if sys_beg_gene is None or sys_end_gene is None:
            return self._invalid_coord(
                cluster_id,
                [],
                f"Missing '{col_map['start_gene']}' or '{col_map['end_gene']}' columns in systems TSV",
            )

        try:
            gene_list = [
                g.strip()
                for g in str(getattr(row, col_map["genes_in_cluster"])).split(",")
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
        verbose: bool = False,
    ) -> typing.Iterator[Cluster]:
        """Extract nucleotide sequences for gene clusters.

        Streams FASTA file to minimize memory usage.

        Args:
            output: Output sink for writing records.
            verbose: Enable detailed logging.

        Returns:
            List of (cluster_id, sequence_length, fasta_file) tuples.
        """
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

        if verbose:
            self.console.print(
                f"[bold blue]{'Processing':>12}[/] {len(valid_coords)} clusters across {num_contigs} contigs/seqs for [bold cyan]{str(self.context.genome_id)}[/]"
            )

        # results = []

        for record in read_fasta(
            self.context.genome_fasta, tar_cache=self._tar_cache
        ):
            if record.id not in contig_groups:
                continue

            if verbose:
                self.console.print(
                    f"[bold blue]{'Loading':>12}[/] contig [blue]{record.id}[/] ({len(contig_groups[record.id])} clusters) for [bold cyan]{self.context.genome_id}[/]"
                )

            for coord in contig_groups[record.id]:
                subseq = str(record.seq[coord.start_coord - 1 : coord.end_coord])
                yield Cluster(coord.cluster_id, subseq, source=coord.fasta_file)

                if verbose:
                    self.console.print(
                        f"[bold blue]{'Extracted':>12}[/] [cyan]{coord.cluster_id}[/] ({len(subseq)} bp) for [bold cyan]{self.context.genome_id}[/]"
                    )

            del contig_groups[record.id]

            if not contig_groups:
                break

        if contig_groups:
            for seq_id in contig_groups:
                self.console.print(
                    f"[bold red]{'Error':>12}[/] Contig {seq_id} not found in genome [bold cyan]{self.context.genome_id}[/]"
                )

    def extract_proteins_from_coordinates(
        self,
        coordinates: typing.List[SystemCoordinates],
        verbose: bool = False,
    ) -> typing.Iterable[Protein]:
        """Extract protein sequences from gene coordinates.

        Args:
            coordinates: List of cluster coordinates.
            output_file: Output file handle for writing sequences.
            verbose: Enable detailed logging.

        Returns:
            Dictionary mapping protein_id to sequence length.
        """
        valid_coords = [c for c in coordinates if c.valid]

        if not valid_coords:
            self.console.print(
                f"[bold yellow]{'Warning':>12}[/] No valid clusters for [bold cyan]{self.context.genome_id}[/]"
            )
            return

        all_gene_ids = set()
        for coord in valid_coords:
            all_gene_ids.update(coord.genes)

        self.protein_idx.load_subset(all_gene_ids)

        total_genes = len(all_gene_ids)

        if verbose:
            self.console.print(
                f"[bold blue]{'Processing':>12}[/] {total_genes} proteins from {len(valid_coords)} clusters for [bold cyan]{str(self.context.genome_id)}[/]"
            )

        n_extracted = 0
        for coord in valid_coords:
            for gene_id in coord.genes:
                if seq := self.protein_idx.get_with_fallback(gene_id):
                    protein_id = f"{coord.cluster_id}__{gene_id}"
                    yield Protein(
                        protein_id,
                        seq,
                        cluster_id=coord.cluster_id,
                    )
                    n_extracted += 1
                else:
                    if verbose:
                        self.console.print(
                            f"[bold yellow]{'Warning':>12}[/] Protein {gene_id} "
                            f"not found for cluster [bold cyan]{coord.cluster_id}[/]"
                        )

        self.console.print(
            f"[bold blue]{'Extracted':>12}[/] {n_extracted} proteins from "
            f"{len(coordinates)} clusters in [bold cyan]{self.context.genome_id}[/]"
        )


class FastaGFFDataset(BaseDataset):
    """Dataset for extracting sequences from FASTA/GFF files."""

    def __init__(
        self,
        inputs: typing.List[pathlib.Path],
        column_mapping: typing.Optional[typing.Dict[str, str]] = None,
    ) -> None:
        """Initialize the FastaGFFDataset class.

        Args:
            inputs: List of input paths (metadata TSV or individual
                files).
            column_mapping: Custom column mapping. If None, uses
                default generic mapping.
        """
        super().__init__()
        self.inputs = inputs

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
        """Create GenomeContext from metadata row.

        Args:
            row: Dictionary with file paths from metadata TSV.
            genome_id: Genome identifier.

        Returns:
            GenomeContext instance.
        """
        return GenomeContext(
            genome_id=genome_id,
            cluster_tsv=pathlib.Path(row.cluster_tsv),
            gff_file=pathlib.Path(row.gff_file),
            genome_fasta=pathlib.Path(row.genome_fasta_file),
            protein_fasta=pathlib.Path(row.protein_fasta_file),
            column_mapping=self.column_mapping,
            system_loader=self._load_and_filter_systems,
        )

    def _load_and_filter_systems(
        self, tsv_path: pathlib.Path, console: Console
    ) -> pandas.DataFrame:
        """Load systems TSV file.

        Args:
            tsv_path: Path to systems TSV.
            console: Rich console for logging.

        Returns:
            Polars DataFrame with system data.
        """
        df = pandas.read_csv(tsv_path, sep="\t")
        return df

    def extract_clusters(
        self,
        progress: rich.progress.Progress,
    ) -> typing.Iterable[Cluster]:
        progress.console.print(
            f"[bold blue]{'Using':>12}[/] cluster metadata file: [magenta]{self.cluster_metadata}[/]"
        )

        df = pandas.read_csv(self.cluster_metadata, sep="\t")

        results = []

        genome_count = 0
        task = progress.add_task(
            f"[bold blue]{'Processing':>9}[/] gene clusters", total=len(df)
        )

        for row in df.itertuples():
            genome_id = getattr(row, "genome_id", f"genome_{genome_count:07}")
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

            with GenomeResources(context, progress.console) as resources:
                coordinates = resources.coordinates
                valid_count = sum(1 for c in coordinates if c.valid)

                progress.update(
                    task,
                    description=f"[bold blue]{'Processing':>9}[/] genome: [bold cyan]{genome_id}[/] - validated {valid_count}/{len(coordinates)} clusters",
                )

                n_extracted = 0
                for cluster in resources.extract_genome_sequences(self.verbose):
                    yield cluster
                    n_extracted += 1

                contig_groups = collections.defaultdict(list)
                for coord in coordinates:
                    if coord.valid:
                        contig_groups[coord.seq_id].append(coord)
                num_contigs = len(contig_groups)

                progress.console.print(
                    f"[bold blue]{'Extracted':>12}[/] {n_extracted} "
                    f"gene clusters across {num_contigs} contigs/seqs for "
                    f"[bold cyan]{context.genome_id}[/]"
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

            progress.update(task, advance=1)

        progress.remove_task(task)
        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(results)} gene clusters "
            f"from {len(df)} genomes/MAGs"
        )

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        cluster_ids: typing.Container[str],
    ) -> typing.Dict[str, int]:

        df = pandas.read_csv(self.cluster_metadata, sep="\t")

        genome_count = 0
        task = progress.add_task(
            f"[bold blue]{'Processing':>9}[/] protein sequences", total=len(df)
        )

        for row in df.itertuples():
            genome_id = getattr(row, "genome_id", f"genome_{genome_count:07}")
            genome_count += 1
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

            with GenomeResources(context, progress.console) as resources:
                coordinates = resources.coordinates
                if cluster_ids:
                    coordinates = [
                        c
                        for c in coordinates
                        if c.cluster_id in cluster_ids
                    ]
                yield from resources.extract_proteins_from_coordinates(
                    coordinates,
                    verbose=self.verbose,
                )

            progress.update(task, advance=1)

        progress.remove_task(task)

    def _log_protein_summary(
        self,
        progress: rich.progress.Progress,
        results: typing.Dict,
        representatives: typing.Optional[typing.Container[str]],
    ):
        """Log summary of protein extraction results.

        Args:
            progress: Rich progress bar instance.
            results: Dictionary of extracted proteins.
            representatives: Container of representative cluster IDs.
        """
        rep_count = "all"
        if representatives:
            if hasattr(representatives, "__len__"):
                rep_count = str(len(representatives))
            else:
                rep_count = "specified"

        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(results):,} proteins from {rep_count} representative gene clusters"
        )
