"""Non-overwriting future artifact writers; never merge numerical methods."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import uuid
import numpy as np
from .spectrum import FullRecordFFTPSD

class FFTArtifactError(RuntimeError):
    pass

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def write_fft_product(*, result: FullRecordFFTPSD, output_directory: Path, channel: str, input_provenance: dict) -> tuple[Path, Path]:
    """Future-only numerical product writer; refuses overwrite and identifies FFT method."""
    output_directory.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    base = output_directory / f"fft-full-record-{channel}-{token}"
    npz, record = base.with_suffix(".npz"), base.with_suffix(".json")
    if npz.exists() or record.exists():
        raise FFTArtifactError("refusing to overwrite FFT artifact")
    np.savez_compressed(npz, frequency_hz=result.frequency_hz, psd_v2_per_hz=result.psd_v2_per_hz, asd_v_per_sqrt_hz=result.asd_v_per_sqrt_hz)
    metadata = {"method_identity": "fft_full_record", "channel": channel, "unit": "V", "psd_unit": "V^2/Hz", "asd_unit": "V/sqrt(Hz)", "artifact_sha256": _sha256(npz), "input_provenance": input_provenance}
    record.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return npz, record

def combined_display_specification(*, fft_artifact: Path, c_lpsd_artifact: Path) -> dict:
    """Return display-only provenance; never returns a merged numerical spectrum."""
    return {"kind": "display_only_combination", "fft_method_identity": "fft_full_record", "c_method_identity": "c_lpsd", "fft_artifact": str(fft_artifact), "c_lpsd_artifact": str(c_lpsd_artifact), "numerical_products_merged": False}
