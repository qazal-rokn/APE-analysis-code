# Provisional one-command real APE analysis (29 October 2026)

This integration **calls the existing frozen functions without modifying scientific algorithms**.

Run from the repository checkout on the user's Mac (Python environment containing APE dependencies and the original vendor package):

```bash
git switch development/likelihood
git pull --ff-only
export PYTHONPATH="$PWD/src:$PWD"
python -B -m tools.ape_real_one_command --config configs/ape_real_one_command.example.json --dry-run
python -B -m tools.ape_real_one_command --config configs/ape_real_one_command.example.json
```

If the native C stage has finished and only downstream inference stopped, rerun the **same configuration** with `--resume`; it verifies input and PS hashes and does not recalculate native C. An interrupted native-C run without a completed checkpoint requires a separate fresh run directory, to avoid reusing unverified partial C output.

Use an existing frozen calibration `stage_result.json` or configure `source.run_config` to execute the existing macOS DXD frozen calibration first. The latter is a frozen-segment method, not general arbitrary voltage-series cleaning. No manual intermediate file copying, manual C output joining, or likelihood execution is required.

Default sample config runs the **original verified C executable** in the existing `lpsd-linux` Docker container. It verifies the frozen HDF5 stage result and binary SHA, copies calibrated HDF5 only when missing and verifies its hash, uses existing C argument assembly to run FFT and Constant-Q batches for both independent quadratures, streams the C *native PS* column (rad²), and then calls the existing vendor fitter, vendor MLE, local-candidate code, and asymptotic expected 95% upper-limit function.

For matching **already existing native PS files**, use `native_ps.mode="reuse"`, provide `calibrated_hdf5_sha256` equal to the frozen calibration hash, and provide both quadratures under `spectra` with absolute path, SHA-256, `unit="rad^2"` and `schema="frequency_PS_2_columns"`. This disallows treating historical PSD (rad²/Hz) as PS. It also explicitly prevents reusing older spectra with a different calibration hash.

The resulting `downstream/frequency_results.csv`, `downstream/summary.json`, and `downstream/expected_coupling_diagnostic.png` are provisional. Raw and frozen data and upstream algorithms stay untouched. There is **no additional f/Qa or N^(-1/4) factor**. C averaging and vendor noise distribution are not replaced.

Important limitations:
- The sample's `quadrature_power_fraction=[0.5,0.5]` is provisional; it is **not experimentally approved physical QWP X/Y normalization**. All optical-response inputs must be reviewed before physical interpretation.
- The available integrated vendor MLE bridge is **single-bin**. There is no scientifically validated real APE multibin halo/C-window inference stage; the script cannot silently manufacture one.
- The product is an **expected** limit and local candidate ranking, **not** observed 95% exclusion, globally calibrated significance or publication-grade sensitivity.
- The sample path settings refer to the previously recorded Mac/Docker environment and may not exist in a new environment.
- The `--dry-run` checks the selected source/binary but does **not** perform a full statistical numerical regression.
- An unrestricted five-million-frequency MLE scan can take extensive CPU time. No speed promise or validated duration.
- If a stage fails or an output already exists, it preserves the failure state rather than overwriting prior results. The downstream chain supports checksum-verified restart using `--resume`. Native C does not resume partially completed batches automatically.

This adds the missing software orchestration to the extent supported by *existing* methods. It does NOT complete the unfinished physical multibin or observed-limit methodology and must not be called a validated physical exclusion pipeline.
