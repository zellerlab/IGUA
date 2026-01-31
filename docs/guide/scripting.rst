Scripting IGUA
==============

To facilitate using IGUA on arbitrary datasets, IGUA since ``v0.2.0`` provides
an API which can be used to easily pass various configurations to IGUA.


Basic usage
-----------

.. code:: python

    import tempfile
    import rich.progress

    from igua.pipeline import ClusteringPipeline
    from igua.dataset.genbank import GenBankDataset

    dataset = GenBankDataset("clusters.gbk")
    with tempfile.TemporaryDirectory() as tempdir:
        with rich.progress.Progress() as progress:
            pipeline = ClusteringPipeline(tempdir, jobs=8, progress=progress)
            results = pipeline.run(dataset=dataset)
            

Combining datasets
------------------

.. code:: python

    from igua.dataset.genbank import GenBankDataset

    dataset = GenBankDataset("clusters.gbk")


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