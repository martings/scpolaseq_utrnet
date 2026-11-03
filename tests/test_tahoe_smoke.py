import json
from pathlib import Path

from scpolaseq_utrnet.tahoe_smoke import run_fixture


def test_checked_in_fixture_runs_end_to_end():
    root = Path(__file__).resolve().parents[1]
    fixture = root / "examples" / "tahoe_smoke_fixture.json"
    result = run_fixture(fixture)
    assert result["status"] == "pass"
    assert result["cells"] == 2
    assert result["raw_count_sum_after_marker_exclusion"] == 8
    assert result["biological_claim"] is False
