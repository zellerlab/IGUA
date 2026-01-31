API Reference
=============

This section contains a comprehensive reference of the public API of IGUA.

.. currentmodule:: igua

Pipeline
--------

.. currentmodule:: igua.pipeline

.. autosummary::

   ClusteringParameters
   ClusteringPipeline
   ClusteringResult

.. toctree::
   :caption: Data Structures
   :maxdepth: 1
   :hidden:

   Pipeline <pipeline>


Clustering Strategy
-------------------

.. currentmodule:: igua.clustering

.. autosummary::

   ClusteringStrategy
   HierarchicalClustering
   LinearClustering

.. toctree::
   :caption: Data Structures
   :maxdepth: 1
   :hidden:

   Clustering <clustering>


Datasets
--------

.. toctree::
   :caption: Datasets
   :maxdepth: 1
   :hidden:

   Base <dataset/base>
   List <dataset/list>
   GenBank <dataset/genbank>
   antiSMASH <dataset/antismash>
   FASTA/GFF <dataset/fasta_gff>
   DefenseFinder <dataset/defensefinder>

Base
^^^^

.. currentmodule:: igua.dataset.base

.. autosummary::

    BaseDataset


List
^^^^

.. currentmodule:: igua.dataset.list

.. autosummary::

    DatasetList


GenBank
^^^^^^^

.. currentmodule:: igua.dataset.genbank

.. autosummary::

    GenBankDataset


antiSMASH
^^^^^^^^^

.. currentmodule:: igua.dataset.antismash

.. autosummary::

    AntiSMASHGenBankDataset
    AntiSMASHZipDataset


Fasta/GFF
^^^^^^^^^

.. currentmodule:: igua.dataset.fasta_gff

.. autosummary::

    FastaGFFDataset


DefenseFinder
^^^^^^^^^^^^^

.. currentmodule:: igua.dataset.defensefinder

.. autosummary::

    DefenseFinderDataset