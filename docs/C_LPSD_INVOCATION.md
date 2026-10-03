# C-LPSD invocation

Every C call uses an argument vector, an explicit temporary `LPSDCFN` containing
`FORMAT fDNRb`, `-D /strain`, `-h 0`, `-w -2`, `-p 200`, `-r 1`, `-l 50`,
`-a 100`, `-m 1`, `-x 1`, and a configured memory power. `tmax=(N-1)/fs` is
explicit because the C HDF5 defaults path does not infer dataset length.

Synthetic benchmarks use deliberately small `Jdes` values only to measure and
validate behavior. They are not production-grid claims. Persistent artifacts are
created under `output_directory/synthetic-benchmarks/run-<uuid>` and never
overwrite an earlier run.
