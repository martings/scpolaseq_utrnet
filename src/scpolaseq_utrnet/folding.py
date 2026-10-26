"""Optional cached RNAplfold unpaired probability backend."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import pandas as pd


def annotate_accessibility(regulatory_sites: pd.DataFrame, isoforms: pd.DataFrame,
                           cache: Path, *, window: int = 80, span: int = 40,
                           unpaired: int = 5) -> pd.DataFrame:
    result = regulatory_sites.copy()
    if result.empty:
        result["accessibility"] = pd.Series(dtype=float)
        result["accessibility_status"] = pd.Series(dtype=str)
        return result
    if not (0 < unpaired <= span <= window):
        raise ValueError("RNAplfold requires 0 < unpaired <= span <= window")
    executable = shutil.which("RNAplfold")
    result["accessibility"] = float("nan")
    result["accessibility_status"] = "missing_backend" if executable is None else "unavailable"
    if executable is None:
        return result
    cache.mkdir(parents=True, exist_ok=True)
    for iso in isoforms.itertuples(index=False):
        if not iso.sequence:
            continue
        key = hashlib.sha256(json.dumps([iso.sequence, window, span, unpaired]).encode()).hexdigest()
        target = cache / f"{key}.json"
        if target.exists():
            probabilities = json.loads(target.read_text())
        else:
            with tempfile.TemporaryDirectory() as td:
                proc = subprocess.run([executable, "-W", str(window), "-L", str(span), "-u", str(unpaired)],
                                      input=f">sequence\n{iso.sequence}\n", text=True,
                                      cwd=td, capture_output=True, timeout=180, check=False)
                if proc.returncode:
                    continue
                lunp = Path(td) / "sequence_lunp"
                if not lunp.exists():
                    continue
                probabilities = {}
                for line in lunp.read_text().splitlines():
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.split()
                    try:
                        probabilities[int(fields[0])] = float(fields[unpaired])
                    except (IndexError, ValueError):
                        continue
                target.write_text(json.dumps(probabilities))
        mask = result.isoform_id == iso.isoform_id
        for idx in result.index[mask]:
            start = int(result.at[idx, "isoform_start"])
            end = int(result.at[idx, "isoform_end"])
            positions = range(start + 1, end + 1)
            values = [probabilities.get(str(p), probabilities.get(p)) for p in positions]
            if any(v is None for v in values):
                result.at[idx, "accessibility_status"] = "boundary_or_missing_probability"
            else:
                result.at[idx, "accessibility"] = sum(values) / len(values)
                result.at[idx, "accessibility_status"] = "predicted_RNAplfold"
    return result


def common_region_accessibility(regulatory_sites: pd.DataFrame) -> pd.DataFrame:
    """Compare the same genomic binding region across reconstructed isoforms."""
    columns = ["regulator", "chrom", "genomic_start", "genomic_end", "strand",
               "n_isoforms", "min_accessibility", "max_accessibility", "range_accessibility"]
    if regulatory_sites.empty or "accessibility" not in regulatory_sites:
        return pd.DataFrame(columns=columns)
    observed = regulatory_sites[regulatory_sites.accessibility.notna()]
    if observed.empty:
        return pd.DataFrame(columns=columns)
    key = ["regulator", "chrom", "genomic_start", "genomic_end", "strand"]
    result = observed.groupby(key, as_index=False).agg(
        n_isoforms=("isoform_id", "nunique"), min_accessibility=("accessibility", "min"),
        max_accessibility=("accessibility", "max"))
    result = result[result.n_isoforms > 1].copy()
    result["range_accessibility"] = result.max_accessibility - result.min_accessibility
    return result[columns]
