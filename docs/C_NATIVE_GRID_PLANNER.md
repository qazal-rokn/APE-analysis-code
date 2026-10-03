# C-Native Grid Planner

This approved no-data planner selects the densest *unchanged C method-0* grid inside the fixed 65--9000 Hz search interval for the validated segment (`N=308,717,000`, `fs=20,000 Hz`). It changes no C estimator setting.

It mirrors `lpsd.c` exactly for the C command-line grid:

```text
f_j = fmin exp[j log(fmax/fmin)/(Jdes-1)]
L_j = round(fs / { f_j [exp(log(fmax/fmin)/(Jdes-1)) - 1] })
```

The planner chooses the largest integer `Jdes` whose `L_j` values all fit the segment, then validates every node. The longest requirement is at the first node. It writes an immutable NPZ containing `frequency_hz` and `required_dft_samples`, with an adjacent JSON run record containing all selected C options and a comparison to the archived APE grid. It never reads PEM values, writes HDF5, or invokes `lpsd-exec`.

Run:

```sh
.venv/bin/python tools/plan_c_native_grid.py --help
.venv/bin/python tools/plan_c_native_grid.py --config configs/run.yaml
```

The manifest records a direct exact C invocation (`-n Jdes -N 0`). The C executable uses `start = nspec * iter`; therefore a fixed `-n` batching plan is exact only if that `n` divides `Jdes`. This planner does not silently change a batch size or use partial final batches.

The planner additionally validates the C-reported Kaiser `SBIN=8.15` absolute lower bound using the unchanged C check equation. It is provenance and a refusal guard only; the approved search range remains 65--9000 Hz.
