"""Metric registry compatibility gates; no historical metric values are recomputed."""
from __future__ import annotations

import pytest

from scpolaseq_utrnet.metric_registry import (
    contrast_delta,
    contrast_record,
    interpretation_allowed,
    legacy_metric_definition,
    metric_definition,
    order_pas_by_proximity,
)


def test_two_pas_order_on_positive_and_negative_strands():
    plus = [{"site_id": "p-distal", "start": 200, "strand": "+"},
            {"site_id": "p-proximal", "start": 100, "strand": "+"}]
    minus = [{"site_id": "m-distal", "start": 100, "strand": "-"},
             {"site_id": "m-proximal", "start": 200, "strand": "-"}]
    assert order_pas_by_proximity(plus) == ["p-proximal", "p-distal"]
    assert order_pas_by_proximity(minus) == ["m-proximal", "m-distal"]


def test_three_or_more_pas_preserve_strand_aware_rank_order():
    plus = [{"site_id": f"p{i}", "start": pos, "strand": "+"}
            for i, pos in enumerate([300, 100, 200])]
    minus = [{"site_id": f"m{i}", "start": pos, "strand": "-"}
             for i, pos in enumerate([100, 300, 200])]
    assert order_pas_by_proximity(plus) == ["p1", "p2", "p0"]
    assert order_pas_by_proximity(minus) == ["m1", "m2", "m0"]


def test_explicit_contrast_orientation_handles_positive_and_negative_delta():
    assert contrast_delta("pas_usage_proportion", 0.2, 0.8,
                         orientation="group_b_minus_group_a") == pytest.approx(0.6)
    assert contrast_delta("pas_usage_proportion", 0.8, 0.2,
                         orientation="group_b_minus_group_a") == pytest.approx(-0.6)
    assert contrast_delta("multisite_shortening_index", 0.2, 0.8,
                         orientation="group_a_minus_group_b") == pytest.approx(-0.6)
    record = contrast_record(
        "pas_usage_proportion", "treated", 0.2, "control", 0.8,
        orientation="group_a_minus_group_b", source_implementation="fixture-v1",
    )
    assert record["delta"] == pytest.approx(-0.6)
    assert record["group_a"]["label"] == "treated"
    assert record["orientation"] == "group_a_minus_group_b"


def test_new_contrast_requires_supported_explicit_orientation():
    with pytest.raises(TypeError, match="orientation"):
        contrast_delta("pas_usage_proportion", 0.2, 0.8)
    with pytest.raises(ValueError, match="orientation must"):
        contrast_delta("pas_usage_proportion", 0.2, 0.8, orientation="legacy_unspecified")


def test_legacy_metric_values_are_gated_without_renaming_or_rewriting():
    entry = legacy_metric_definition("si")
    assert entry is not None
    assert entry["metric_id"] == "multisite_shortening_index"
    assert entry["source_field"] == "si"
    assert entry["compatibility_status"] == "definition_verified"
    assert interpretation_allowed(entry["metric_id"])
    assert legacy_metric_definition("unknown_apa_field") is None
    fixed = metric_definition("fixed_basis_genomic_magnetization")
    assert fixed["orientation"] == "distal_high_proximal_low"
    assert fixed["range"] == [-1.0, 1.0]
    proximal = metric_definition("proximal_pas_fraction")
    assert proximal["orientation"] == "proximal_high_distal_low; 1=proximal"
    assert not interpretation_allowed("differential_apa_statistic")


def test_pas_order_requires_explicit_single_strand_and_missing_site_fields():
    with pytest.raises(ValueError, match="one explicit strand"):
        order_pas_by_proximity([{"site_id": "a", "start": 1, "strand": "+"},
                                {"site_id": "b", "start": 2, "strand": "-"}])
    with pytest.raises(ValueError, match="site_id and start"):
        order_pas_by_proximity([{"site_id": "a", "strand": "+"}])
