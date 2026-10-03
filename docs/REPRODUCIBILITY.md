# Reproducibility

The project-local `.venv` was created with Python 3.12. Installed primary
packages: NumPy 2.5.3, h5py 3.16.0, PyYAML 6.0.3, Matplotlib 3.11.2.
Run `.venv/bin/python -m pip freeze` to capture the full environment.

The C executable path is configuration data. Runs must record its SHA-256,
resolved path, exact argument vector, controlled `LPSDCFN`, stdout/stderr, C
configuration, dynamic-library details, and output hashes.

