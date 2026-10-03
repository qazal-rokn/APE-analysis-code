# FFT Extension Architecture

`src/ape_c_lpsd/fft_extension/` is intentionally isolated from both the
validated segment extractor and the C bridge:

- `config.py`: YAML contract and range/window/detrend validation only.
- `preflight.py`: no-data full-record memory and range refusal.
- `spectrum.py`: future one-channel full-record one-sided FFT kernel.
- `artifacts.py`: non-overwriting future method-specific products and
  display-only combination provenance.

The interface accepts `samples_v`, `sample_rate_hz`, and validated metadata
from an upstream owner. It cannot locate, extract, or clean APE source data.
