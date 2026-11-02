"""Generic regulatory-signature/v0.1 fixtures and missingness tests."""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from scpolaseq_utrnet.signature import validate_regulatory_signature, write_regulatory_signature


ROOT = Path(__file__).resolve().parents[1]


def _fixture(name: str):
    return json.loads((ROOT / "examples" / name).read_text(encoding="utf-8"))


def test_pbmc_and_nonimmune_signatures_validate_without_hardcoded_context():
    pbmc = _fixture("regulatory_signature_pbmc.json")
    adipocyte = _fixture("regulatory_signature_adipocyte.json")
    assert validate_regulatory_signature(pbmc) == []
    assert validate_regulatory_signature(adipocyte) == []
    assert adipocyte["dataset"]["species"] == "Mus musculus"
    assert "timepoint" in adipocyte["observations"][0]
    assert pbmc["modalities"]["small_rna"]["status"] == "not_measured"


def test_signature_round_trip_is_canonical(tmp_path: Path):
    payload = _fixture("regulatory_signature_adipocyte.json")
    shuffled = dict(payload)
    shuffled["measurements"] = list(reversed(payload["measurements"]))
    first = write_regulatory_signature(payload, tmp_path / "first.json")
    second = write_regulatory_signature(shuffled, tmp_path / "second.json")
    assert first.read_bytes() == second.read_bytes()
    assert validate_regulatory_signature(json.loads(first.read_text())) == []


def test_measured_zero_is_distinct_from_missing_zero_coverage():
    payload = _fixture("regulatory_signature_adipocyte.json")
    measured_zero = next(item for item in payload["measurements"] if item["status"] == "measured")
    missing = next(item for item in payload["measurements"] if item["status"] == "missing_zero_coverage")
    assert measured_zero["value"] == 0
    assert measured_zero["status"] == "measured"
    assert missing["value"] is None and missing["coverage"] == 0

    missing["value"] = 0
    assert any("null value" in error for error in validate_regulatory_signature(payload))


def test_signature_rejects_unknown_observation_and_unprovenanced_measurement():
    payload = _fixture("regulatory_signature_pbmc.json")
    payload["measurements"][0]["observation_id"] = "not-in-context"
    payload["measurements"][0].pop("source")
    errors = validate_regulatory_signature(payload)
    assert any("unknown observation_id" in error for error in errors)
    assert any("requires source provenance" in error for error in errors)


@pytest.mark.parametrize("mutation", ["value", "modalities", "context"])
def test_runtime_validation_enforces_published_schema(mutation):
    payload = _fixture("regulatory_signature_pbmc.json")
    if mutation == "value":
        payload["measurements"][0].pop("value")
    elif mutation == "modalities":
        payload.pop("modalities")
    else:
        payload["observations"][0].pop("context")
    errors = validate_regulatory_signature(payload)
    assert any("JSON Schema" in error for error in errors)


@pytest.mark.parametrize("nonfinite", [math.nan, math.inf, -math.inf])
def test_nonfinite_nested_values_are_rejected_without_touching_existing_file(tmp_path: Path, nonfinite):
    payload = _fixture("regulatory_signature_adipocyte.json")
    payload["measurements"][0]["statistical_support"] = {"nested": [nonfinite]}
    destination = tmp_path / "signature.json"
    destination.write_text("previous complete artifact")
    with pytest.raises(ValueError, match="non-finite"):
        write_regulatory_signature(payload, destination)
    assert destination.read_text() == "previous complete artifact"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["signature.json"]


def test_packaged_and_published_json_schemas_are_identical():
    published = json.loads((ROOT / "docs/contracts/regulatory-signature-v0.1.schema.json").read_text())
    packaged = json.loads((ROOT / "src/scpolaseq_utrnet/schemas/regulatory-signature-v0.1.schema.json").read_text())
    assert packaged == published
