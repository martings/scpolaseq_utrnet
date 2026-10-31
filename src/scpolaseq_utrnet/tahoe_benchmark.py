"""Transparent zero-delta and regularized linear baselines for Tahoe."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


def zero_delta(n_samples: int, n_genes: int) -> np.ndarray:
    if n_samples < 0 or n_genes < 0:
        raise ValueError("output dimensions must be nonnegative")
    return np.zeros((n_samples, n_genes), dtype=np.float32)


@dataclass
class RidgeDeltaPredictor:
    """Multi-output ridge regression over caller-encoded basal/context/dose features."""

    alpha: float = 1.0
    coefficients_: np.ndarray | None = None
    intercept_: np.ndarray | None = None

    def fit(self, features: np.ndarray, targets: np.ndarray) -> "RidgeDeltaPredictor":
        x = np.asarray(features, dtype=np.float64)
        y = np.asarray(targets, dtype=np.float64)
        if x.ndim != 2 or y.ndim != 2 or x.shape[0] != y.shape[0] or x.shape[0] < 2:
            raise ValueError("features and targets must be aligned 2D arrays with at least two rows")
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("training arrays must contain only finite values")
        if self.alpha <= 0:
            raise ValueError("alpha must be positive")
        self._x_mean = x.mean(axis=0)
        self.intercept_ = y.mean(axis=0)
        centered_x = x - self._x_mean
        centered_y = y - self.intercept_
        gram = centered_x @ centered_x.T
        dual = np.linalg.solve(gram + self.alpha * np.eye(len(x)), centered_y)
        self.coefficients_ = centered_x.T @ dual
        return self

    def predict(self, features: np.ndarray) -> np.ndarray:
        if self.coefficients_ is None or self.intercept_ is None:
            raise RuntimeError("fit must be called before predict")
        x = np.asarray(features, dtype=np.float64)
        if x.ndim != 2 or x.shape[1] != self.coefficients_.shape[0]:
            raise ValueError("features have the wrong shape")
        return (x - self._x_mean) @ self.coefficients_ + self.intercept_


def held_out_condition_mask(conditions: list[str], held_out: set[str]) -> np.ndarray:
    if not held_out:
        raise ValueError("at least one complete condition must be held out")
    condition_array = np.asarray(conditions, dtype=str)
    if not held_out.issubset(set(condition_array)):
        raise ValueError("held-out condition absent from input")
    return np.isin(condition_array, sorted(held_out))


def response_metrics(
    truth: np.ndarray,
    prediction: np.ndarray,
    *,
    top_k: int = 50,
) -> dict[str, Any]:
    """Report change error, sign agreement and top-|delta| gene recovery."""
    y = np.asarray(truth, dtype=np.float64)
    p = np.asarray(prediction, dtype=np.float64)
    if y.shape != p.shape or y.ndim != 2 or not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError("truth and prediction must be aligned finite 2D arrays")
    if top_k < 1:
        raise ValueError("top_k must be positive")
    direction_mask = y != 0
    direction = float(np.mean(np.sign(y[direction_mask]) == np.sign(p[direction_mask]))) if direction_mask.any() else float("nan")
    k = min(top_k, y.shape[1])
    true_top = np.argpartition(np.abs(y), -k, axis=1)[:, -k:]
    pred_top = np.argpartition(np.abs(p), -k, axis=1)[:, -k:]
    recovery = float(np.mean([
        len(set(a.tolist()) & set(b.tolist())) / k for a, b in zip(true_top, pred_top)
    ])) if len(y) else float("nan")
    return {
        "mse": float(np.mean((y - p) ** 2)) if y.size else float("nan"),
        "directional_agreement_nonzero_truth": direction,
        "top_k_de_recovery": recovery,
        "top_k": k,
        "n_samples": int(y.shape[0]),
        "n_genes": int(y.shape[1]),
    }


def per_condition_metrics(
    truth: np.ndarray,
    prediction: np.ndarray,
    conditions: list[str],
    *,
    top_k: int = 50,
) -> pd.DataFrame:
    y = np.asarray(truth)
    p = np.asarray(prediction)
    labels = np.asarray(conditions, dtype=str)
    if len(labels) != y.shape[0]:
        raise ValueError("one condition label is required per sample")
    records = []
    for label in sorted(set(labels)):
        mask = labels == label
        records.append({"condition": label, **response_metrics(y[mask], p[mask], top_k=top_k)})
    return pd.DataFrame.from_records(records)
