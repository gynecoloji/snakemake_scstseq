"""Platform readers on in-memory vendor directories (10x HDF5 written with
h5py by the same writer the executable test case uses). No scanpy needed."""
import importlib.util
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import load_sample  # noqa: E402
import scst_common as sc_  # noqa: E402

_spec = importlib.util.spec_from_file_location("make_testdata", ROOT / ".test" / "make_testdata.py")
mtd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mtd)

GENES = ["GENE_A", "GENE_B", "MT-CO1"]


def _xenium_dir(tmp_path, cells_format="csv"):
    d = tmp_path / "xenium"
    d.mkdir()
    X = np.array([[5, 0, 1, 2, 0], [0, 3, 0, 0, 1], [1, 1, 1, 0, 0]])  # 3 cells x (3 genes + 2 controls)
    names = GENES + ["NegControlProbe_00001", "UnassignedCodeword_0001"]
    types = ["Gene Expression"] * 3 + ["Negative Control Probe", "Unassigned Codeword"]
    mtd.write_10x_h5(d / "cell_feature_matrix.h5", X, ["c1", "c2", "c3"], names, ["id1", "id2", "id3", "nc1", "uc1"], types)
    cells = pd.DataFrame(
        {"cell_id": ["c1", "c2", "c3"], "x_centroid": [10.0, 20.0, 30.0], "y_centroid": [1.0, 2.0, 3.0],
         "transcript_counts": [8, 4, 3], "cell_area": [50.0, 60.0, 70.0], "nucleus_area": [15.0, 18.0, 21.0]}
    )
    if cells_format == "csv":
        cells.to_csv(d / "cells.csv.gz", index=False)
    else:
        cells.to_parquet(d / "cells.parquet", index=False)
    return d, X


def test_read_10x_h5_orientation(tmp_path):
    d, X = _xenium_dir(tmp_path)
    a = sc_.read_10x_h5(d / "cell_feature_matrix.h5")
    assert a.shape == (3, 5)
    assert a.X.toarray().tolist() == X.tolist()
    assert list(a.obs_names) == ["c1", "c2", "c3"]
    assert a.var["feature_types"].tolist()[:3] == ["Gene Expression"] * 3


def test_read_xenium_counts_and_removes_controls(tmp_path):
    d, X = _xenium_dir(tmp_path)
    a = load_sample.read_xenium(str(d))
    assert list(a.var_names) == GENES
    assert a.obs["control_counts"].tolist() == [2, 1, 0]
    assert a.obs["control_frac"].tolist() == pytest.approx([2 / 8, 1 / 4, 0.0])
    assert a.obsm["spatial"].tolist() == [[10.0, 1.0], [20.0, 2.0], [30.0, 3.0]]
    assert a.obs["cell_area"].tolist() == [50.0, 60.0, 70.0]
    assert a.X.toarray().tolist() == X[:, :3].tolist()
    assert "Negative Control Probe" in a.uns["scst"]["control_feature_types"]


def test_read_xenium_parquet(tmp_path):
    pytest.importorskip("pyarrow")
    d, _ = _xenium_dir(tmp_path, cells_format="parquet")
    a = load_sample.read_xenium(str(d))
    assert a.obsm["spatial"][1].tolist() == [20.0, 2.0]


def test_read_xenium_missing_cell_raises(tmp_path):
    d, _ = _xenium_dir(tmp_path)
    df = pd.read_csv(d / "cells.csv.gz")
    df[df["cell_id"] != "c2"].to_csv(d / "cells.csv.gz", index=False)
    with pytest.raises(ValueError, match="absent from the cells table"):
        load_sample.read_xenium(str(d))


def _visium_dir(tmp_path, header=True):
    d = tmp_path / "visium"
    (d / "spatial").mkdir(parents=True)
    X = np.array([[5, 0, 1], [0, 3, 0], [1, 1, 1], [2, 2, 2]])
    barcodes = ["AAAA-1", "CCCC-1", "GGGG-1", "TTTT-1"]
    mtd.write_10x_h5(d / "filtered_feature_bc_matrix.h5", X, barcodes, GENES, ["i1", "i2", "i3"], ["Gene Expression"] * 3, library_id="v")
    rows = [
        ["AAAA-1", 1, 0, 0, 100, 200],
        ["CCCC-1", 1, 0, 2, 100, 300],
        ["GGGG-1", 1, 1, 1, 187, 250],
        ["TTTT-1", 0, 5, 5, 900, 900],   # in the matrix but out of tissue
        ["NNNN-1", 0, 6, 6, 999, 999],   # not in the matrix
    ]
    cols = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if header:
        pd.DataFrame(rows, columns=cols).to_csv(d / "spatial" / "tissue_positions.csv", index=False)
    else:
        pd.DataFrame(rows).to_csv(d / "spatial" / "tissue_positions_list.csv", index=False, header=False)
    (d / "spatial" / "scalefactors_json.json").write_text(json.dumps({"spot_diameter_fullres": 45.0, "tissue_hires_scalef": 0.1}))
    return d


@pytest.mark.parametrize("header", [True, False])
def test_read_visium(tmp_path, header):
    d = _visium_dir(tmp_path, header=header)
    a = load_sample.read_visium(str(d), "v")
    assert list(a.obs_names) == ["AAAA-1", "CCCC-1", "GGGG-1"]  # out-of-tissue spot dropped
    assert a.uns["scst"]["n_out_of_tissue"] == 1
    # (x, y) = (pxl_col, pxl_row)
    assert a.obsm["spatial"].tolist() == [[200.0, 100.0], [300.0, 100.0], [250.0, 187.0]]
    assert a.obs["array_row"].tolist() == [0, 0, 1]
    assert a.uns["spatial"]["v"]["scalefactors"]["spot_diameter_fullres"] == 45.0
    assert a.uns["spatial"]["v"]["images"] == {}


def test_read_visium_missing_barcode_raises(tmp_path):
    d = _visium_dir(tmp_path)
    pos = pd.read_csv(d / "spatial" / "tissue_positions.csv")
    pos[pos["barcode"] != "CCCC-1"].to_csv(d / "spatial" / "tissue_positions.csv", index=False)
    with pytest.raises(ValueError, match="absent from the tissue positions"):
        load_sample.read_visium(str(d), "v")


def test_load_sample_prefixes_barcodes(tmp_path):
    d, _ = _xenium_dir(tmp_path)
    a = load_sample.load_sample("xenium", str(d), "s1", "tumor")
    assert list(a.obs_names) == ["s1_c1", "s1_c2", "s1_c3"]
    assert a.obs["barcode"].tolist() == ["c1", "c2", "c3"]
    assert a.obs["condition"].unique().tolist() == ["tumor"]
    assert a.uns["scst"]["platform"] == "xenium" and a.uns["scst"]["coord_unit"] == "µm"


def test_load_sample_unknown_platform(tmp_path):
    with pytest.raises(ValueError):
        load_sample.load_sample("cosmx", str(tmp_path), "s", "")


def test_make_testdata_h5_roundtrip(tmp_path):
    """The fixture generator's HDF5 writer and the workflow reader agree."""
    X = np.array([[0, 0], [1, 2]])
    mtd.write_10x_h5(tmp_path / "m.h5", X, ["b1", "b2"], ["g1", "g2"], ["i1", "i2"], ["Gene Expression"] * 2)
    a = sc_.read_10x_h5(tmp_path / "m.h5")
    assert a.X.toarray().tolist() == X.tolist()
