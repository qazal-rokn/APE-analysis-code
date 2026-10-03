# C-LPSD production readiness and batching report

Date: 2026-09-26. This is a read-only estimate based only on the completed
calibrated pilot, frozen grid manifest `c-native-grid-1a37f2795f5f415a84656d06b4769e15`, and existing file metadata. It does not authorize production.

## Frozen inputs

- Segment: 308,717,000 samples at 20,000 Hz, with the frozen three-part automatic manifest.
- Grid: 4,947,016 native method-0 nodes, 65--9000 Hz.
- C profile: `-h 0 -w -2 -p 200 -r 1 -l 50 -a 100 -m 1 -x 1 -M 20`, output `fDNRb`.
- Calibrated inputs: separate dimensionless X and Y ellipticity HDF5 datasets.

## Completed pilot measurements

| Region | Node batch | X runtime (s) | Y runtime (s) | Peak child RSS (KB) | Text output per channel |
|---|---:|---:|---:|---:|---:|
| Low | 0 | 280.397 | 281.061 | 214,712,320 | 2,325 B |
| Middle | 309188 | 54.044 | 53.867 | 214,712,320 | 2,348 B |
| High | 618376 | 33.248 | 33.266 | 214,712,320 | 2,356 B |

Each batch contained eight exact frozen-grid nodes. Per-node runtimes are therefore X/Y low 35.05/35.13 s, middle 6.76/6.73 s, and high 4.16/4.16 s. Peak RSS was reported cumulatively by the platform child-process accounting; it is an observed upper bound, not a clean per-job memory measurement.

## Conservative runtime estimate

A linear extrapolation is uncertain: C LPSD cost is strongly frequency-dependent and pilot batches are small. For both channels and all 9,894,032 node evaluations:

- optimistic high-region rate: about 41.2 million s (477 days);
- middle-region rate: about 66.7 million s (772 days);
- simple equal-third low/middle/high regional average: about 151.8 million s (1,757 days, 4.8 years);
- pessimistic low-region rate: about 347.3 million s (4,020 days).

These are *single-worker-equivalent estimates*, not a production forecast. They exclude queueing, I/O contention, retries, checksums, and post-processing. The low-frequency region dominates risk. More approved pilot profiling would be required before treating any estimate as operationally reliable.

## Storage estimate

The two validated input HDF5 files occupy 2,472,020,816 B each (4.94 GB total). The eight-node pilot output is approximately 2.34 kB per channel, or about 293 B/node. At that observed textual `fDNRb` density, production text output is about 1.45 GB/channel, 2.90 GB for X+Y, before logs, parsed arrays, checksums, and filesystem overhead. Reserve at least 10 GB for inputs, C text outputs, logs, and non-overwriting checkpoint artifacts; this is an estimate, not a free-space preflight approval.

## Exact batching recommendation

`Jdes = 4,947,016 = 8 * 618,377`. The only useful exact fixed batch plan is **618,377 nodes per batch, 8 batches per channel** (`-n 618377`, `-N 0` through `-N 7`). It produces 16 independent X/Y jobs. The C source starts a batch at `nspec * iter`; because 618,377 divides Jdes, the jobs are non-overlapping and cover the same frozen global nodes exactly once. Changing only `-n` and `-N` under this divisor rule changes execution partitioning, not frequency nodes, calibration, or scientific settings.

The eight jobs are not balanced in runtime: low-frequency batches are expected to be much slower. A scheduler should allow independent job resources and conservative wall-times rather than assuming equal completion.

## Future checkpoint/resume layout

A future approved production directory must be unique, e.g. `results/c-lpsd-production-<uuid>/`, and contain immutable `manifest/`, `inputs/`, `commands/`, `logs/`, `raw-output/x/`, `raw-output/y/`, and `parsed/` subdirectories. Each completed batch writes one output, stdout, stderr, command JSON, parsed data, and SHA-256. Resume must scan only for a fully validated five-column `fDNRb` output with the exact expected index range and text-precision-aware node identity; otherwise it must rerun only that missing/invalid batch into a new attempt path. Never append output files or overwrite a completed artifact.

## Exact future commands (not executed)

For channel `x_ellipticity`, dataset `/x_ellipticity`, batch `B` in `0..7`:

```text
LPSDCFN=<production>/commands/x-B.cfg lpsd-exec -d -D /x_ellipticity -f 20000 -h 0 -i <pilot-or-approved-production-input>/x_ellipticity.h5 -o <production>/raw-output/x/batch-B.txt -b 0 -e 15435.84995 -s 65 -t 9000 -J 4947016 -n 618377 -N B -w -2 -p 200 -r 1 -l 50 -a 100 -m 1 -x 1 -M 20
```

For `y_ellipticity`, replace dataset, input, and output path with `/y_ellipticity`, `y_ellipticity.h5`, and `raw-output/y/batch-B.txt`. Resume uses the identical command only for batches failing the saved strict output/provenance validation.

No production command was launched while preparing this report.
