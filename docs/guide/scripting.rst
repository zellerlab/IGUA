Scripting IGUA
==============

To facilitate using IGUA on arbitrary datasets, IGUA since ``v0.2.0`` provides
an API which can be used to easily pass various configurations to IGUA, 
including complex or custom datasets.


Basic usage
-----------

The `~igua.pipeline.Pipeline` class implements the same features as the 
command line. It can be used to cluster the gene clusters provided in 
a dataset. The following snippet shows the minimum code needed to cluster 
gene clusters contained in a GenBank file:

.. code:: python

    import rich.progress

    from igua.pipeline import Pipeline
    from igua.dataset.genbank import GenBankDataset

    dataset = GenBankDataset("clusters.gbk")
    with rich.progress.Progress() as progress:
        pipeline = Pipeline(jobs=8, progress=progress)
        results = pipeline.run(dataset)
        print(results.gcfs)
            

Combining datasets
------------------

IGUA provides multiple dataset classes that can be used to manage IGUA
inputs. 

.. code:: python

    from igua.dataset.genbank import GenBankDataset

    dataset = GenBankDataset("clusters.gbk")


The `~igua.dataset.list.DatasetList` class can be used to combine several
datasets inside a single dataset class. For instance, to combine gene 
clusters from two different GenBank files, the following code can be 
used:

.. code:: python

    from igua.dataset.genbank import GenBankDataset
    from igua.dataset.list import DatasetList

    dataset = DatasetList([
        GenBankDataset("clusters1.gbk"),
        GenBankDataset("clusters2.gbk"),
    ])


Implementing a custom dataset class
-----------------------------------

IGUA tries its best to support various formats, but there are always edge 
cases. What if your data is in a ``tar`` archive? What if you don't have 
gene translations inside a GenBank file but in a supplementary FASTA file?
You can implement your own dataset! 

A `~igua.dataset.base.BaseDataset` needs to implement only two methods:

- `~igua.dataset.base.BaseDataset.extract_clusters` to extract the gene 
  genomic sequences of the gene clusters to process.
- `~igua.dataset.base.BaseDataset.extract_proteins` to extract the protein
  sequences of the gene clusters to process.


For instance, let's implement a `BaseDataset` using the output tables produced
by `GECCO <https://gecco.embl.de>`_: not the GenBank files, but the 
``clusters.tsv`` and ``genes.tsv`` table.

Initialization
^^^^^^^^^^^^^^

Let's start by defining a new subclass, and adding the paths we need to 
the constructor. We can parse the tables directly in the ``__init__`` 
method:

.. code:: python

    import pandas
    from pathlib import Path
    from igua.data.base import BaseDataset

    class GECCOClustersDataset(BaseDataset):

        def __init__(self, genome_fna, clusters_tsv, genes_tsv):
            self.genome_fna = genome_fna
            self.clusters_tsv = clusters_tsv
            self.genes_tsv = genes_tsv
            self.clusters = pandas.read_table(self.clusters_tsv, sep="\t")
            self.genes = pandas.read_table(self.genes_tsv, sep="\t")


Cluster extraction
^^^^^^^^^^^^^^^^^^

Now we need to implement the ``extract_clusters`` method, which yields the 
clusters to process with their genomic sequence. To do so, we can use Biopython 
to stream records from the FASTA file, and use the clusters TSV

.. code:: python

    import Bio.SeqIO
    from igua.data.base import BaseDataset, Cluster

    class GECCOClustersDataset(BaseDataset):

        # ... __init__ ...

        def extract_clusters(self, progress):
            for record in Bio.SeqIO.parse(self.genome_fna, "fasta"):
                # extract the clusters predicted in this record and yield
                # a Cluster object for each of them
                record_clusters = self.clusters_tsv[self.clusters_tsv["sequence_id"] == record.id]
                for row in record_clusters.itertuples():
                    yield Cluster(
                        id=row.cluster_id,
                        sequence=str(record.seq[row.start-1:row.end]),
                        source=str(self.genome_fna)
                    )

Note that we are not ensuring that all clusters have been extracted successfully, 
which may be a desirable check to add to the code. You could use a dictionary
to track which clusters in ``self.clusters_tsv['cluster_id']`` have been 
found. 

Protein extraction
^^^^^^^^^^^^^^^^^^

Now moving on to the ``extract_proteins`` method. This method takes an 
additional argument, ``clusters``, which is a container (typically a list)
of cluster IDs from which to extract proteins from. 


.. code:: python

    import Bio.SeqIO
    from igua.data.base import BaseDataset, Protein

    class GECCOClustersDataset(BaseDataset):

        # ... __init__ ...

        # ... extract_clusters ...

        def extract_proteins(self, progress, clusters):
            # subset the clusters table with only the IDs of the clusters to
            # extract, which are the representatives of the previous stages
            clusters_to_extract = self.clusters_tsv[self.clusters_tsv["cluster_id"].str.isin(clusters)]
            for record in Bio.SeqIO.parse(self.genome_fna, "fasta"):
                # extract the clusters predicted in this record and extract
                # the genes corresponding to this cluster
                record_clusters = clusters_to_extract[clusters_to_extract["sequence_id"] == record.id]
                for cluster_row in record_clusters.itertuples():
                    cluster_gene_ids = cluster_row.proteins.split(";")
                    cluster_genes = self.genes_tsv[self.genes_tsv["protein_id"].str.isin(cluster_genes_ids)]
                    # yield a protein for every gene, using Biopython to 
                    # translate the genomic sequence to get the protein
                    for gene_row in cluster_genes.itertuples():
                        seq = str(record.seq[gene_row.start-1:gene_row.end].translate())
                        yield Protein(
                            id=gene_row.protein_id,
                            sequence=seq,
                            cluster_id=cluster_row.cluster_id,
                        )


With this method implemented, we have a working dataset which can load data 
described by a GECCO table. It's not the best in terms of performance 
(we could for instance index ``self.genes`` by ``protein_id`` to accelerate
the retrieval of table rows), but should be enough to get us started!

Using the dataset
^^^^^^^^^^^^^^^^^

Now that we have a working dataset class, we can simply instantiate the
dataset and pass it to a pipeline:

.. code:: python

    dataset = GECCOClustersDataset("genome.fna", "clusters.tsv", "genes.tsv")
    with rich.progress.Progress() as progress:
        pipeline = Pipeline(jobs=8, progress=progress)
        results = pipeline.run(dataset)