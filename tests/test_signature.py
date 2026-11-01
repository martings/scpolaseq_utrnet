"""Generic regulatory-signature/v0.1 fixtures and missingness tests."""
from __future__ import annotations

import json
from pathlib import Path

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
