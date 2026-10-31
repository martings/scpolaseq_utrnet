import numpy as np
import pytest

from scpolaseq_utrnet.tahoe_benchmark import (
    RidgeDeltaPredictor,
    held_out_condition_mask,
    response_metrics,
    cross_validate_ridge,
)


def test_zero_delta_and_ridge_baselines_and_response_metrics():
    x = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    y = np.column_stack([2 * x[:, 0] - 1, x[:, 1] - 0.5])
    model = RidgeDeltaPredictor(alpha=0.01).fit(x, y)
    pred = model.predict(x)
    assert pred.shape == y.shape
    assert response_metrics(y, pred, top_k=1)["mse"] < 0.01


def test_condition_holdout_uses_whole_conditions():
    mask = held_out_condition_mask(["a", "b", "a", "c"], {"b"})
    assert mask.tolist() == [False, True, False, False]
    with pytest.raises(ValueError, match="absent"):
        held_out_condition_mask(["a", "b"], {"missing"})


def test_zero_and_ridge_share_complete_condition_holdout():
    x = np.array([[0.0], [0.1], [1.0], [1.1]])
    y = np.array([[0.0, 0.0], [0.1, 0.0], [1.0, 0.0], [1.1, 0.0]])
    result = cross_validate_ridge(x, y, ["train", "train", "held", "held"], {"held"}, top_k=1)
    assert result["train_rows"] == 2
    assert result["test_rows"] == 2
    assert set(result["per_condition"]) == {"zero_delta", "ridge"}
