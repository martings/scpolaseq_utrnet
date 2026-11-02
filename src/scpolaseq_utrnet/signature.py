"""Small generic `regulatory-signature/v0.1` data contract and canonical writer."""
from __future__ import annotations

from copy import deepcopy
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any

from jsonschema import Draft202012Validator

SIGNATURE_VERSION = "regulatory-signature/v0.1"
FEATURE_CLASSES = {
    "gene_expression", "differential_expression", "pas_usage", "apa_metric",
    "transcript_association", "regulatory_region_availability",
    "regulatory_association", "functional_enrichment", "pathway_activity",
    "perturbation_response",
}
MEASUREMENT_STATUSES = {
    "measured", "experimental_observation", "annotation_derived", "candidate_inferred",
    "missing_unmeasured", "missing_zero_coverage", "ambiguous",
}
MODALITY_STATUSES = {"measured", "partially_measured", "not_measured"}
_SCHEMA_PATH = Path(__file__).parent / "schemas" / "regulatory-signature-v0.1.schema.json"


def _nonfinite_paths(value: Any, path: str = "$") -> list[str]:
    if isinstance(value, float) and not math.isfinite(value):
        return [path]
    if isinstance(value, dict):
        return [found for key, item in value.items() for found in _nonfinite_paths(item, f"{path}.{key}")]
    if isinstance(value, list):
        return [found for index, item in enumerate(value) for found in _nonfinite_paths(item, f"{path}[{index}]")]
    return []


def validate_regulatory_signature(payload: dict[str, Any]) -> list[str]:
    """Validate semantic distinctions and flexible context without disease assumptions."""
    errors: list[str] = []
    if not isinstance(payload, dict):
        return ["signature payload must be an object"]
    if not _SCHEMA_PATH.is_file():
        return [f"published signature schema is unavailable: {_SCHEMA_PATH}"]
    try:
        schema = json.loads(_SCHEMA_PATH.read_text(encoding="utf-8"))
        Draft202012Validator.check_schema(schema)
        validator = Draft202012Validator(schema)
        errors.extend(
            f"JSON Schema {list(error.absolute_path)}: {error.message}"
            for error in validator.iter_errors(payload)
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        errors.append(f"cannot load published signature schema: {exc}")
    errors.extend(f"non-finite number at {path}" for path in _nonfinite_paths(payload))
    if payload.get("schema_version") != SIGNATURE_VERSION:
        errors.append(f"schema_version must be {SIGNATURE_VERSION!r}")
    if not isinstance(payload.get("dataset"), dict) or not payload["dataset"].get("dataset_id"):
        errors.append("dataset.dataset_id is required")
    observations = payload.get("observations")
    if not isinstance(observations, list):
        errors.append("observations must be a list")
        observations = []
    observation_ids = [str(item.get("observation_id", "")) for item in observations if isinstance(item, dict)]
    if len(observation_ids) != len(observations) or any(not item for item in observation_ids):
        errors.append("each observation requires observation_id")
    if len(observation_ids) != len(set(observation_ids)):
        errors.append("observation_id values must be unique")
    for observation in observations:
        if not isinstance(observation, dict):
            continue
        if not isinstance(observation.get("context", {}), dict):
            errors.append(f"observation {observation.get('observation_id')} context must be an object")

    modalities = payload.get("modalities", {})
    if not isinstance(modalities, dict):
        errors.append("modalities must be an object keyed by modality")
        modalities = {}
    for name, modality in modalities.items():
        if (not isinstance(modality, dict) or not isinstance(modality.get("status"), str)
                or modality.get("status") not in MODALITY_STATUSES):
            errors.append(f"modality {name!r} needs a valid status")
            continue
        if modality["status"] == "not_measured" and modality.get("measurement_scale") not in (None, ""):
            errors.append(f"unmeasured modality {name!r} cannot claim a measurement scale")

    measurements = payload.get("measurements", [])
    if not isinstance(measurements, list):
        errors.append("measurements must be a list")
        measurements = []
    measurement_ids: set[str] = set()
    for index, item in enumerate(measurements):
        if not isinstance(item, dict):
            errors.append(f"measurement {index} must be an object")
            continue
        mid = str(item.get("measurement_id", ""))
        if not mid:
            errors.append(f"measurement {index} requires measurement_id")
        elif mid in measurement_ids:
            errors.append(f"duplicate measurement_id {mid!r}")
        measurement_ids.add(mid)
        feature_class = item.get("feature_class")
        if not isinstance(feature_class, str) or feature_class not in FEATURE_CLASSES:
            errors.append(f"measurement {mid!r} has unknown feature_class")
        if not item.get("feature_id") or not item.get("entity_id"):
            errors.append(f"measurement {mid!r} requires feature_id and entity_id")
        status = item.get("status")
        if not isinstance(status, str) or status not in MEASUREMENT_STATUSES:
            errors.append(f"measurement {mid!r} has unknown status")
        if isinstance(status, str) and status in {"missing_unmeasured", "missing_zero_coverage"} and item.get("value") is not None:
            errors.append(f"measurement {mid!r} must use null value for explicit missing status")
        if status == "missing_zero_coverage" and item.get("coverage") != 0:
            errors.append(f"measurement {mid!r} missing_zero_coverage requires coverage=0")
        if status == "missing_unmeasured" and item.get("coverage") == 0:
            errors.append(f"measurement {mid!r} has zero coverage; use missing_zero_coverage")
        observation_id = item.get("observation_id")
        if observation_id is not None and str(observation_id) not in set(observation_ids):
            errors.append(f"measurement {mid!r} references unknown observation_id")
        if not item.get("source") or not isinstance(item.get("source"), dict):
            errors.append(f"measurement {mid!r} requires source provenance")
        if feature_class in {"functional_enrichment", "pathway_activity"}:
            if not item.get("program_id") or not item.get("program_source") or not item.get("program_version"):
                errors.append(f"functional measurement {mid!r} requires program_id, program_source, and program_version")
            if feature_class == "functional_enrichment" and "gene_universe_n" not in item:
                errors.append(f"functional enrichment {mid!r} requires gene_universe_n")
    return errors


def write_regulatory_signature(payload: dict[str, Any], path: str | Path) -> Path:
    """Validate and write a stable JSON signature with deterministic measurement order."""
    errors = validate_regulatory_signature(payload)
    if errors:
        raise ValueError("Invalid regulatory-signature/v0.1:\n- " + "\n- ".join(errors))
    canonical = deepcopy(payload)
    canonical["measurements"] = sorted(
        canonical.get("measurements", []),
        key=lambda item: (
            str(item.get("feature_class", "")), str(item.get("feature_id", "")),
            str(item.get("observation_id", "")), str(item.get("measurement_id", "")),
        ),
    )
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(canonical, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    temporary: str | None = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=destination.parent,
                                        prefix=f".{destination.name}.", suffix=".tmp", delete=False) as handle:
            temporary = handle.name
            handle.write(serialized)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)
        raise
    return destination
