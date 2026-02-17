Running IGUA
============

IGUA is a method for clustering gene clusters. It assumes that you have 
both the genomic sequences for each gene cluster, as well as gene calling
annotations with the translation of each protein. 

Basic usage (GenBank input)
---------------------------

The most common format to pass gene cluster records to IGUA is the GenBank
format. In GenBank inputs, each complete record is taken as a gene cluster,
and protein sequences are extracted from the ``/translation`` qualifiers
of ``CDS`` features.

To run IGUA on a GenBank file named ``clusters.gbk`` and output the 
resulting gene cluster families in a table named ``gcfs.tsv``, run in 
the console:

.. code:: console

    $ igua -i clusters.gbk -o gcfs.tsv

.. caution::

    The GenBank locus identifier will be used as the name of each gene cluster. 
    This may cause problems with gene clusters obtained with some tools, 
    such as antiSMASH. If the input contains duplicate identifiers, the first 
    gene cluster with a given identifier will be used, and a warning will 
    be displayed.

IGUA can take more than one input file, and will just concatenate
records from different files if given more than one input:

.. code:: console

    $ igua -i clusters1.gbk -i clusters2.gbk ...

.. caution::

    Beware of the order of the clusters! Because of the heuristics used in 
    MMSeqs2, the order of the gene clusters in the inputs actually influences
    the GCFs returned by IGUA. You may get different GCFs for the same gene
    clusters ordered differently! Make sure the gene clusters are always 
    given in a consistent order to avoid non-deterministic results.


Understanding the output
------------------------

The main output of IGUA is a TSV file which assigns a Gene Cluster Family to 
each gene cluster found in the input. The GCF identifiers are arbitrary, and
the prefix can be changed with the ``--prefix`` flag. The table also records
the original file from which each record was obtained to facilitate resource
management. The table is written to the filename given with the ``--output`` 
flag. The output table contains the following columns:

``cluster_id``
    The identifier of the gene cluster.

``gcf_id``
    The GCF this gene cluster belongs to, assigned to during the third stage
    of the pipeline (protein composition clustering step).

``nucleotide_representative``
    The representative gene cluster this gene cluster was assigned to during the
    second stage of the pipeline (nucleotide clustering step).

``fragment_representative``
    The representative gene cluster this gene cluster was assigned to during the
    first stage of the pipeline (nucleotide deduplication step).

``filename``
    The filename this gene cluster originates from. This can be useful to trace
    the source of a gene cluster which may have been obtained from different 
    input files.


The sequences of the representative proteins extracted from each cluster 
can be saved to a FASTA file with the ``--features`` flag. These proteins are
used for compositional representation of gene clusters, and can be used to
transfer annotations to the GCF representatives. The final compositional matrix 
for each GCF representative, which can be useful for computing distances 
between GCFs, can be saved as an `anndata` HDF5 matrix to a filename given 
with the `--compositions` flag.


Resource management
-------------------

MMseqs needs a fast scratch space to work with intermediate files while running
linear clustering. By default, this will use a temporary folder obtained with
``tempfile.TemporaryDirectory``, which typically lies inside ``/tmp``. To use a 
different folder, use the ``--workdir`` flag. 

The temporary folder will at some point of the run contain a local copy of 
the gene cluster sequences in FASTA 2-line format for use by MMSeqs2, so it
needs to be large enough to accomodate all sequences. This may be a problem
for some systems where ``/tmp`` is generally small. If you run out of 
storage space, consider setting ``--workdir`` to another location with more
storage available.


Tuning IGUA
-----------

IGUA can be configured in many ways that are not shown by default. To see 
all parameters of IGUA, run it with the ``--help-all`` flag:

.. code:: console

    $ igua --help-all


Nucleotide deduplication and clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first two stages of IGUA rely on using MMSeqs with parameters that have been
tuned during IGUA's development. Nevertheless, they can be modified with 
the various ``--dedup-*`` and ``--nuc-*`` flags of the command line if more
control is desired.


Protein composition clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the protein composition clustering step of IGUA will use 
hierarchical clustering with **average** linkage clustering and a relative
distance threshold of :math:`0.8`, which corresponds to clusters indide a 
GCF having at most 20% of esimtated difference at the amino-acid level. These
two options can be changed with the ``--clustering-method`` and 
``--clustering-distance`` flags. 

The pairwise distance matrix is computed using a weighted Manhattan distance
for sparse matrices. By default, the weight of each protein corresponds to its
length in amino-acids, which means larger proteins influence the distance 
more. This may not be desired in cases where small proteins may have a higher
functional importance (e.g. precursor peptides in RiPP BGCs), and can be 
turned of with the ``--no-clustering-weight`` option if desired.

Additionally, the precision of the pairwise distance matrix used for the clustering
can be lowered to reduce memory usage, using ``single`` or ``half`` precision
float point numbers instead of the ``double`` precision used by default. Use the 
``--clustering-precision`` flag to control numerical precision. *Note that 
half-precision may lead to very degraded accuracy of the clustering*.



