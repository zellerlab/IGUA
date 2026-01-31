Scripting IGUA
==============

To facilitate using IGUA on arbitrary datasets, IGUA since ``v0.2.0`` provides
an API which can be used to easily pass various configurations to IGUA.


Basic usage
-----------

The `~igua.pipeline.Pipeline` class implements the same features as the 
command line. It can be used to cluster the gene clusters provided in 
a dataset. The following snippet shows how to cluster gene clusters 
in a GenBank file:

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


