# Configuration B: EXPERIMENTAL halo / multi-bin diagnostics ONLY

Configuration B is **not** part of the Configuration A statistical analysis.
It contains exploratory astrophysical/line-shape, Kaiser-window, adjacent-bin
correlation, and provisional QWP-power template studies.

Existing diagnostic source code remains in `../ape_post_lpsd/`, including:
- `ape_axion_halo_kernel_diagnostic.py`
- `ape_native_kaiser_kernel.py`
- `ape_provisional_qwp_multibin_vendor.py`
- `ape_signal_weighted_covariance_diagnostic.py`
- `physics/axion_halo.py`

This folder does not run automatically and intentionally provides **no**
production Configuration B command. A multi-bin physical likelihood requires
a checked estimator covariance/coverage model and a validated signal response.
Outputs from these scripts must never be relabeled as Configuration A results.

**Rule:** Configuration A stays on the vendor single-bin likelihood path;
Configuration B experiments never modify its files or inputs.
