API Reference
=============

This section contains a comprehensive reference of the public API of IGUA.

.. currentmodule:: igua

Pipeline
--------

.. currentmodule:: igua.pipeline

.. autosummary::
   :nosignatures:

   PipelineParameters
   PipelineResult
   Pipeline

.. toctree::
   :caption: Pipeline
   :maxdepth: 1
   :hidden:

   Pipeline <pipeline>


Clustering Strategy
-------------------

.. currentmodule:: igua.clustering

.. autosummary::
   :nosignatures:

   ClusteringStrategy
   HierarchicalClustering
   LinearClustering

.. toctree::
   :caption: Clustering Strategy
   :maxdepth: 1
   :hidden:

   Clustering <clustering>


MMSeqs Driver
-------------

.. currentmodule:: igua.mmseqs

.. autosummary::
   :nosignatures:

   MMSeqs

.. toctree::
   :caption: MMSeqs Driver
   :maxdepth: 1
   :hidden:

   MMSeqs <mmseqs>



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
   :nosignatures:

    Cluster
    Protein
    BaseDataset


List
^^^^

.. currentmodule:: igua.dataset.list

.. autosummary::
   :nosignatures:

    DatasetList


GenBank
^^^^^^^

.. currentmodule:: igua.dataset.genbank

.. autosummary::
   :nosignatures:

    GenBankDataset


antiSMASH
^^^^^^^^^

.. currentmodule:: igua.dataset.antismash

.. autosummary::
   :nosignatures:

    AntiSMASHGenBankDataset
    AntiSMASHZipDataset


Fasta/GFF
^^^^^^^^^

.. currentmodule:: igua.dataset.fasta_gff

.. autosummary::
   :nosignatures:

    FastaGFFDataset


DefenseFinder
^^^^^^^^^^^^^

.. currentmodule:: igua.dataset.defensefinder

.. autosummary::
   :nosignatures:

    DefenseFinderDataset