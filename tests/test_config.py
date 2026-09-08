"""Validate config/config.yaml against the JSON schema and sanity-check samples.csv."""
import csv
import pathlib

import pytest
import yaml
from jsonschema import Draft7Validator

ROOT = pathlib.Path(__file__).resolve().parent.parent

VALID_PLATFORMS = {"xenium", "visium"}


def load_schema():
    return yaml.safe_load((ROOT / "workflow/schemas/config.schema.yaml").read_text())


def check_samples_sheet(path):
    rows = list(csv.DictReader(path.open()))
    assert rows, f"{path} has no data rows"
    for col in ("sample_id", "platform", "path"):
        assert col in rows[0], f"{path} is missing the {col} column"
    ids = [r["sample_id"].strip() for r in rows]
    assert all(ids), f"empty sample_id found in {path}"
    assert len(ids) == len(set(ids)), f"duplicate sample_id found in {path}"
    for r in rows:
        assert r["platform"] in VALID_PLATFORMS, (
            f"{path}: sample {r['sample_id']} has platform {r['platform']!r}, "
            f"expected one of {sorted(VALID_PLATFORMS)}"
        )
        assert r["path"].strip(), f"{path}: sample {r['sample_id']} has an empty path"


def test_schema_is_valid():
    Draft7Validator.check_schema(load_schema())


def test_config_matches_schema():
    config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    Draft7Validator(load_schema()).validate(config)


def test_every_schema_property_has_a_default_or_is_an_object():
    """The catalog renders defaults from the schema; scalar leaves need one."""
    def walk(node, path):
        for name, sub in node.get("properties", {}).items():
            here = f"{path}.{name}" if path else name
            if sub.get("type") == "object":
                walk(sub, here)
            else:
                assert "default" in sub, f"schema property {here} has no default"

    walk(load_schema(), "")


def test_samples_sheet():
    check_samples_sheet(ROOT / "config/samples.csv")


def test_test_config_matches_schema():
    test_config = ROOT / ".test/config/config.yaml"
    if not test_config.exists():
        pytest.skip(".test/ catalog case not present")
    Draft7Validator(load_schema()).validate(yaml.safe_load(test_config.read_text()))


def test_test_samples_sheet():
    test_samples = ROOT / ".test/config/samples.csv"
    if not test_samples.exists():
        pytest.skip(".test/ catalog case not present")
    check_samples_sheet(test_samples)
