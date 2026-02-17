import pandas as pd
import anndata
import pytest


test_cases = ["ecoli", "ecoli_gzip"]


@pytest.mark.parametrize("ref_data_id", test_cases)
def test_regression_gcfs(ref_data_id):
    """Compare output against known good reference"""
    # using defaults (see makefile):
    gcfs = pd.read_csv(f"tests/fixtures/{ref_data_id}_gcfs.tsv", sep="\t")
    gcfs_alternative = pd.read_csv(
        f"tests/fixtures/{ref_data_id}_gcfs_alternative.tsv", sep="\t"
    )
    current_gcfs = pd.read_csv(f"tests/test_output/{ref_data_id}_gcfs.tsv", sep="\t")

    # same shape and number of unique GCFs
    assert len(gcfs) == len(current_gcfs)
    assert gcfs["gcf_id"].nunique() == current_gcfs["gcf_id"].nunique()

    # clustering should be deterministic (if using same parameters)
    # therefore the number of unique values per column should be the same
    pd.testing.assert_series_equal(current_gcfs.nunique(), gcfs.nunique())
    # and the output gcfs should be identical
    try:
        pd.testing.assert_frame_equal(
            gcfs.sort_values("cluster_id").reset_index(drop=True),
            current_gcfs.sort_values("cluster_id").reset_index(drop=True),
        )
    except AssertionError:
        pd.testing.assert_frame_equal(
            gcfs_alternative.sort_values("cluster_id").reset_index(drop=True),
            current_gcfs.sort_values("cluster_id").reset_index(drop=True),
        )


@pytest.mark.parametrize("ref_data_id", test_cases)
def test_regression_compositions(ref_data_id):
    """Test compositions output structure"""
    comp = anndata.read_h5ad(f"tests/fixtures/{ref_data_id}_compositions.h5ad")
    comp_alternative = anndata.read_h5ad(
        f"tests/fixtures/{ref_data_id}_compositions_alternative.h5ad"
    )
    current_comp = anndata.read_h5ad(
        f"tests/test_output/{ref_data_id}_compositions.h5ad"
    )

    assert (
        comp.shape == current_comp.shape or comp_alternative.shape == current_comp.shape
    )
    assert all(comp.var_names == current_comp.var_names) or all(
        comp_alternative.var_names == current_comp.var_names
    )
    assert all(comp.obs_names == current_comp.obs_names) or all(
        comp_alternative.obs_names == current_comp.obs_names
    )
