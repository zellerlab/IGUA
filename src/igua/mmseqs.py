import enum
import os
import pathlib
import subprocess
import typing

import pandas
import rich.progress
import tempfile


COMMANDS = {
  "easy-search",
  "easy-linsearch",
  "easy-cluster",
  "easy-linclust",
  "easy-taxonomy",
  "easy-rbh",
  "search",
  "linsearch",
  "map",
  "rbh",
  "linclust",
  "cluster",
  "clusterupdate",
  "taxonomy",
  "databases",
  "createdb",
  "createindex",
  "createlinindex",
  "convertmsa",
  "tsv2db",
  "tar2db",
  "msa2profile",
  "compress",
  "decompress",
  "rmdb",
  "mvdb",
  "cpdb",
  "lndb",
  "unpackdb",
  "touchdb",
  "createsubdb",
  "concatdbs",
  "mergedbs",
  "subtractdbs",
  "convertalis",
  "createtsv",
  "convert2fasta",
  "result2flat",
  "createseqfiledb",
  "taxonomyreport",
  "extractorfs",
  "extractframes",
  "orftocontig",
  "reverseseq",
  "translatenucs",
  "translateaa",
  "splitsequence",
  "masksequence",
  "extractalignedregion",
  "swapresults",
  "result2rbh",
  "result2msa",
  "result2dnamsa",
  "result2stats",
  "filterresult",
  "offsetalignment",
  "proteinaln2nucl",
  "result2repseq",
  "sortresult",
  "summarizealis",
  "summarizeresult",
  "createtaxdb",
  "createbintaxonomy",
  "addtaxonomy",
  "taxonomyreport",
  "filtertaxdb",
  "filtertaxseqdb",
  "aggregatetax",
  "aggregatetaxweights",
  "lcaalign",
  "lca",
  "majoritylca",
  "multihitdb",
  "multihitsearch",
  "besthitperset",
  "combinepvalperset",
  "mergeresultsbyset",
  "prefilter",
  "ungappedprefilter",
  "kmermatcher",
  "kmersearch",
  "align",
  "alignall",
  "transitivealign",
  "rescorediagonal",
  "alignbykmer",
  "clust",
  "clusthash",
  "mergeclusters",
  "result2profile",
  "msa2result",
  "msa2profile",
  "profile2pssm",
  "profile2consensus",
  "profile2repseq",
  "convertprofiledb",
  "enrich",
  "result2pp",
  "profile2cs",
  "convertca3m",
  "expandaln",
  "expand2profile",
  "view",
  "apply",
  "filterdb",
  "swapdb",
  "prefixid",
  "suffixid",
  "renamedbkeys",
  "diffseqdbs",
  "summarizetabs",
  "gff2db",
  "maskbygff",
  "convertkb",
  "summarizeheaders",
  "nrtotaxmapping",
  "extractdomains",
  "countkmer",
}


class MMSeqs(object):
    """A wrapper around an ``mmseqs`` binary and common parameters.
    """

    def __init__(
        self,
        binary: str = "mmseqs",
        progress: typing.Optional[rich.progress.Progress] = None,
        threads: typing.Optional[int] = None,
        tempdir: typing.Optional[str] = None,
    ) -> None:
        self.binary = binary
        self.progress = progress
        self.threads = threads
        self.tempdir = tempdir or tempfile.TemporaryDirectory().name

    def _wrap_progress(
        self,
        process: subprocess.Popen,
    ) -> int:
        """Wrap the progress output from ``mmseqs`` into a `rich` progress bar.
        """
        assert self.progress is not None
        assert process.stdout is not None

        buffer = bytearray()
        command = ""
        task = self.progress.add_task(f"[bold blue]{'Running':>9}[/] [purple]{command}[/]", total=65)
        bar_column = next(c for c in self.progress.columns if isinstance(c, rich.progress.BarColumn))

        for x in iter(lambda: process.stdout.read(1), b""):
            buffer.append(ord(x))
            if buffer.startswith(b"["):
                # update progress
                self.progress.update(task_id=task, completed=buffer.count(b'='))
            if buffer.endswith(b"\n"):
                # extract current command being run
                _command = next(iter(buffer.split()), b"").decode()
                if _command in COMMANDS:
                    command = _command
                    bar_column.bar_width = 40 - len(command)
                    self.progress.reset(task_id=task, description=f"[bold blue]{'Running':>9}[/] [purple]{command}[/]")
                # clear current buffer
                buffer.clear()

        self.progress.update(task_id=task, visible=False)
        self.progress.remove_task(task)
        return process.wait()

    def run(
        self,
        command: str,
        *args: typing.Union[str, bytes, os.PathLike],
        **kwargs: object
    ) -> subprocess.CompletedProcess:
        """Run an arbitrary ``mmseqs`` command.
        """
        # build actual command line
        params = [self.binary, command, *args]
        for k, v in kwargs.items():
            dash = "-" if len(k) == 1 else "--"
            flag = "".join([dash, k.replace("_", "-")])
            params.append(flag)
            params.append(repr(v))

        # use fixed number of threads
        if self.threads is not None and command in {"linclust", "search"}:
            params.extend(["--threads", str(self.threads)])

        # start mmseqs subprocess
        process = subprocess.Popen(params, stdout=subprocess.PIPE, bufsize=0)

        # wrap progress if a rich progress bar is available
        if self.progress:
            self._wrap_progress(process)
            stdout = stderr = None
        else:
            stdout, stderr = process.communicate()

        # return a completed process instance for compatibility
        return subprocess.CompletedProcess(
            params,
            process.returncode,
        )

    def version(self) -> str:
        try:
            params = [self.binary, "version"]
            process = subprocess.Popen(params, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            return stdout.decode().strip()
        except OSError as err:
            raise RuntimeError(f"Failed to find MMSeqs2 binary {self.binary!r}") from err


class DatabaseType(enum.IntEnum):
    Protein = 1
    Nucleotide = 2


class _MMSeqsFile(object):

    def __init__(self, mmseqs: MMSeqs, path: pathlib.Path):
        self.mmseqs = mmseqs
        self.path = path

    def remove(self):
        """Remove the MMSeqs2 file.
        """
        self.mmseqs.run("rmdb", self.path)


class Database(_MMSeqsFile):

    @classmethod
    def create(
        cls,
        mmseqs: MMSeqs,
        sequences: pathlib.Path,
        path: typing.Optional[pathlib.Path] = None,
    ) -> "Database":
        """Create a new database using a FASTA-formatted sequence file.
        """
        if path is None:
            path = sequences.with_suffix(".db")
        mmseqs.run(
            "createdb",
            sequences,
            path,
            shuffle=0,
            createdb_mode=1,
            write_lookup=0,
            id_offset=0,
            compressed=0,
        ).check_returncode()
        return cls(mmseqs, path)

    def __init__(self, mmseqs: MMSeqs, path: pathlib.Path):
        """Open a database at the given location.
        """
        super().__init__(mmseqs, path)

    def to_fasta(self, path: pathlib.Path) -> pathlib.Path:
        """Convert the sequence database to a two-line FASTA file.
        """
        self.mmseqs.run("convert2fasta", self.path, path).check_returncode()
        return path

    def cluster(
        self,
        output: pathlib.Path,
        *,
        e_value=0.001,
        sequence_identity=0.85,
        coverage=1,
        cluster_mode=0,
        coverage_mode=0,
        spaced_kmer_mode=0,
    ) -> "Clustering":
        """Run linear clustering on the database.
        """
        # run clustering
        self.mmseqs.run(
            "linclust",
            self.path,
            output,
            self.mmseqs.tempdir,
            e=e_value,
            min_seq_id=sequence_identity,
            c=coverage,
            cluster_mode=cluster_mode,
            cov_mode=coverage_mode,
            spaced_kmer_mode=spaced_kmer_mode,
            remove_tmp_files=1,
        ).check_returncode()
        return Clustering(self.mmseqs, output, self)


class Clustering(_MMSeqsFile):

    def __init__(
        self,
        mmseqs: MMSeqs,
        path: pathlib.Path,
        database: Database
    ):
        super().__init__(mmseqs, path)
        self.database = database

    def to_dataframe(self, columns: typing.Tuple[str, str]) -> pandas.DataFrame:
        """Obtain the clustering results as a table.
        """
        with tempfile.NamedTemporaryFile(
            "w",
            suffix=".tsv",
            dir=self.mmseqs.tempdir
        ) as tmp:
            self.mmseqs.run(
                "createtsv",
                self.database.path,
                self.database.path,
                self.path,
                tmp.name,
            ).check_returncode()
            return pandas.read_csv(
                tmp.name,
                sep="\t",
                header=None,
                names=list(columns),
            )

    def to_subdb(self, path: pathlib.Path) -> Database:
        """Generate a sub-database with the representatives of this clustering.
        """
        self.mmseqs.run(
            "createsubdb",
            self.path,
            self.database.path,
            path,
            subdb_mode=1,
        ).check_returncode()
        return Database(self.mmseqs, path)