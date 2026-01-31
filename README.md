# 🦎 IGUA [![Stars](https://img.shields.io/github/stars/zellerlab/IGUA.svg?style=social&maxAge=3600&label=Star)](https://github.com/zellerlab/IGUA/stargazers)

*Iterative Gene clUster Analysis, a high-throughput method for gene cluster family identification.*

[![Actions](https://img.shields.io/github/actions/workflow/status/zellerlab/IGUA/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/zellerlab/IGUA/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/zellerlab/IGUA?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/zellerlab/IGUA/)
[![PyPI](https://img.shields.io/pypi/v/igua.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/igua)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/igua?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/igua)
[![AUR](https://img.shields.io/aur/version/igua?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/igua)
[![Wheel](https://img.shields.io/pypi/wheel/igua.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/igua/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/igua.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/igua/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/igua.svg?logo=python&style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/igua/#files)
[![License](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/zellerlab/igua/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/igua/)
[![GitHub issues](https://img.shields.io/github/issues/zellerlab/IGUA.svg?style=flat-square&maxAge=600)](https://github.com/zellerlab/IGUA/issues)
[![Docs](https://img.shields.io/readthedocs/igua/latest?style=flat-square&maxAge=600)](https://igua.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/zellerlab/IGUA/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/igua?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/igua)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2025.05.15.654203v1)


## 🗺️ Overview

IGUA is a method for high-throughput content-agnostic identification of
Gene Cluster Families (GCFs) from gene clusters of genomic and metagenomic 
origin. It performs three clustering iterations to perform GCF assignment:

- *Fragment mapping identification*: Reduce the input sequence space by 
  identifying which gene clusters are fragments of each other. 
- *Nucleotide deduplication*: Find similar gene clusters in genomic space,
  using linear clustering with lower sequence identity and coverage.
- *Protein representation*: Compute a numerical representation of gene clusters
  in term of protein composition, using representatives from a protein sequence
  clustering, to identify more distant relatives not captured by the previous
  step.

Compared to similar methods such as [BiG-SLiCE](https://github.com/medema-group/bigslice) 
or [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), IGUA does not use Pfam 
domains to represent gene cluster composition, using instead representatives
from an unsupervised clustering. This allows IGUA to accurately account for
proteins that may not be covered by Pfam, and avoids performing a costly annotation
step. The resulting protein representatives can be later annotated indepently
to transfer annotations to the GCFs.


## 🔧 Installing

### Bioconda

IGUA and all of its dependencies are available via [Bioconda](https://anaconda.org/channels/bioconda/packages/igua/overview) and can be installed using e.g., `conda` or `pixi`:

1. First, [set up Bioconda with Pixi or Conda.](https://bioconda.github.io/)

2. Then, install IGUA using the appropriate method:

With `conda`:
```console
$ conda install igua
```

With `pixi`:
```console
$ pixi add igua
```

### Apptainer, Docker, and Singularity

IGUA (and all of its dependencies) can be run using e.g., Docker, Apptainer, and Singularity, using images available [here](https://quay.io/repository/biocontainers/igua?tab=tags&tag=latest).

An example using Apptainer (using IGUA v0.1.0):
```
apptainer pull docker://quay.io/biocontainers/igua:0.1.0--py39h5b94c0b_0
```

### pip

IGUA can be downloaded directly from PyPI, which hosts pre-compiled 
distributions for Linux, MacOS and Windows. Simply install with `pip`:

```console
$ pip install igua
```

**Note that you will need to install MMseqs2 yourself through other means.**

## 💡 Running

### 📥 Inputs

The gene clusters to pass to IGUA must be in GenBank format, with gene 
annotations inside of `CDS` features. Several GenBank files can be passed
to the same pipeline run.

```console
$ igua -i clusters1.gbk -i clusters2.gbk ...
```

The GenBank locus identifier will be used as the name of each gene cluster. This
may cause problems with gene clusters obtained with some tools, such as antiSMASH.
If the input contains duplicate identifiers, the first gene cluster with a given 
identifier will be used, and a warning will be displayed.


### 📤 Outputs

The main output of IGUA is a TSV file which assigns a Gene Cluster Family to 
each gene cluster found in the input. The GCF identifiers are arbitrary, and
the prefix can be changed with the `--prefix` flag. The table will also record
the original file from which each record was obtained to facilitate resource
management. The table is written to the filename given with the `--output` 
flag.

The sequences of the representative proteins extracted from each cluster 
can be saved to a FASTA file with the `--features` flag. These proteins are
used for compositional representation of gene clusters, and can be used to
transfer annotations to the GCF representatives. The final compositional matrix 
for each GCF representative, which can be useful for computing distances 
between GCFs, can be saved as an `anndata` sparse matrix to a filename given 
with the `--compositions` flag.

### 📝 Workspace

MMseqs needs a fast scratch space to work with intermediate files while running
linear clustering. By default, this will use a temporary folder obtained with
`tempfile.TemporaryDirectory`, which typically lies inside `/tmp`. To use a 
different folder, use the `--workdir` flag.

### 🫧 Clustering

By default, IGUA will use **average** linkage clustering and a relative distance 
threshold of `0.8`, which corresponds to clusters inside a GCF having at most
20% of estimated difference at the amino-acid level. These two options can be
changed with the `--clustering-method` and `--clustering-distance` flags.

Additionally, the precision of the distance matrix used for the clustering can
be lowered to reduce memory usage, using `single` or `half` precision floating
point numbers instead of the `double` precision used by default. Use the
`--clustering-precision` flag to control numerical precision.

### ⚙️ Advanced usage

#### 🛠️ Running IGUA on manual gene clusters

It is also possible to provide manually defined gene clusters through a TSV file with the `--dataset-type manual` option. The TSV file must contain at least the following columns:

- (optional) `genome_id`: unique identifier for each genome/metagenome/strain/sample; by default, the row index will be used (_e.g._ `genome_0000000`)
- `cluster_tsv`: path to a TSV file containing the gene cluster definitions for the genome
- `gff_file`: path to the GFF file corresponding to the genome
- `genome_fasta_file`: path to the genome FASTA file
- `protein_fasta_file`: `path to the protein FASTA file

The cluster TSV files must contain at least the following columns:

- `cluster_id`: unique identifier for each gene cluster
- `start_gene`: identifier of the first gene in the cluster
- `end_gene`: identifier of the last gene in the cluster
- `genes_in_cluster`: comma-separated list of gene identifiers in the cluster

```console
$ igua -i input_metadata.tsv --dataset-type manual
```

Alternatively, if the column names in the cluster TSV files differ from the expected ones, a JSON file can be provided with the `--column-mapping` option to map the expected column names to the actual ones. For example:

```json
{
    "cluster_id": "sys_id",
    "start_gene": "sys_beg",
    "end_gene": "sys_end",
    "genes_in_cluster": "protein_in_syst"
}
```

```console
$ igua -i input_metadata.tsv --dataset-type manual --column-mapping column_mapping.json
```

### 🛡️ Running IGUA on DefenseFinder outputs

IGUA can be run on gene clusters obtained with [DefenseFinder](https://github.com/mdmparis/defense-finder) with a dedicated dataset type `--dataset-type defensefinder`. This will parse the DefenseFinder outputs `***_defense_finder_systems.tsv` to extract the gene clusters and their corresponding protein sequences. 

```console
$ igua -i ../defense_finder_metadata.tsv --dataset-type defense-finder
```

This is equivalent to using the manual dataset type with appropriate column mapping, but more convenient. Hence, `defense_finder_metadata.tsv` must contain the following columns:

- (optional) `genome_id`: unique identifier for each genome/metagenome/strain/sample; by default, the row index will be used (_e.g._ `genome_0000000`)
- `defense_systems_tsv`: path to the `***_defense_finder.tsv` file
- `gff_file`: path to the GFF file corresponding to the genome
- `genome_fasta_file`: path to the genome FASTA file
- `protein_fasta_file`: `path to the protein FASTA file
- (optional) `activity`: filter gene clusters by activity (`defense` or `antidefense`)


IGUA can also be run on a single genome by providing all the required files one by one. 

```console 
$ igua \
$  --dataset-type defense-finder \
$   --defense-systems-tsv /path/to/genome-id_defense_finder_systems.tsv \
$   --gff-file /path/to/genome-id.gff \
$   --genome-fasta-file /path/to/genome-id.fa \
$   --protein-fasta--file /path/to/genome-id.faa
```

#### 🛠️ Running IGUA with custom MMseqs2 parameters

Various MMSeqs2 parameters can be adjusted to control the sensitivity and speed
of the clustering steps. For more information, see `igua --help-all`. Please note 
that the default parameters have been optimized for general usage, and changing 
them may lead to suboptimal results.

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/zellerlab/IGUA/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/zellerlab/IGUA/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/zellerlab/IGUA/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).

*This project was developed by [Martin Larralde](https://github.com/althonos/) 
during his PhD project at the [European Molecular Biology Laboratory](https://www.embl.de/) 
and the [Leiden University Medical Center](https://lumc.nl/en/)
in the [Zeller team](https://github.com/zellerlab).*
