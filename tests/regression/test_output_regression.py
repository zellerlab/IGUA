import pandas as pd
import anndata
from pathlib import Path
import pytest


test_cases = ["ecoli", "ecoli_gzip"]


@pytest.mark.parametrize("ref_data_id", test_cases)
def test_gcfs_output_structure(ref_data_id):
    """Test gcfs output maintains expected structure"""
    gcfs = pd.read_csv(f"tests/test_output/{ref_data_id}_gcfs.tsv", sep="\t")

    # file structure tests
    assert set(gcfs.columns) == {
        "cluster_id",
        "cluster_length",
        "gcf_id",
        "gcf_representative",
        "nucleotide_representative",
        "fragment_representative",
        "filename",
    }

    assert len(gcfs) > 0
    assert gcfs["gcf_id"].nunique() <= len(gcfs)  # gcfs <= inputs
    assert all(gcfs["cluster_length"] > 0)
    assert gcfs["gcf_id"].str.startswith("GCF").all()


@pytest.mark.parametrize("ref_data_id", test_cases)
def test_compositions_output(ref_data_id):
    """Test compositions output structure"""
    comp = anndata.read_h5ad(f"tests/test_output/{ref_data_id}_compositions.h5ad")

    assert comp.X.shape[0] > 0  # type: ignore # has gcfs
    assert comp.X.shape[1] > 0  # type: ignore # has features
    # assert "gcf_id" in comp.obs.columns
    # assert "size" in comp.var.columns


@pytest.mark.parametrize("ref_data_id", test_cases)
def test_features_output(ref_data_id):
    """Test protein features FASTA"""
    features_path = Path(f"tests/test_output/{ref_data_id}_features.fa")
    assert features_path.exists()
    assert features_path.stat().st_size > 0

    # count seqs
    with open(features_path) as f:
        seq_count = sum(1 for line in f if line.startswith(">"))
    assert seq_count > 0
