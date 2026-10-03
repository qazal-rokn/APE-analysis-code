# Canonical APE Raw-Data Contract

The first APE pipeline stage is deliberately narrow:

```text
Dewesoft .dxd
    |
    v
five explicitly named recorded channels
    |
    v
bounded-memory synchronized chunks
```

No calibration, filtering, cleaning, interpolation, resampling, trimming, gap
filling, segment selection, spectral estimation, or X/Y combination is allowed
inside this stage.

## Required channels

The run configuration maps exactly five APE roles to five distinct Dewesoft
channel names:

- `x_1f`: recorded first-harmonic X quadrature;
- `y_1f`: recorded first-harmonic Y quadrature;
- `x_2f`: recorded second-harmonic X quadrature;
- `y_2f`: recorded second-harmonic Y quadrature;
- `cavity_transmission`: recorded cavity-transmission monitor.

The reader never guesses these names.

## Chunked production interface

Large Dewesoft files are not retained as one giant in-memory object.
`iter_ape_dxd_chunks` opens the source read-only and yields one synchronized
`APEDXDChunk` at a time. The default chunk size is 1,000,000 samples.

For every chunk the reader verifies:

1. all five arrays and timestamp arrays have the requested length;
2. all selected samples and timestamps are finite;
3. the five channels have identical timestamps;
4. timestamp spacing agrees with Dewesoft's declared sample rate;
5. the first timestamp of each new chunk is continuous with the preceding
   chunk.

No failed check is repaired silently.

The command-line ingestion wrapper uses `inspect_dxd`, which walks the entire
file through this chunked path and retains only the final metadata summary.
The full record is therefore validated without retaining all five channels in
RAM.

`read_ape_dxd` remains only as a compatibility helper for small files where a
fully materialized `APEDXDData` object is genuinely wanted. Large production
APE records should use the chunk iterator.

## Timing rule

For a regular Dewesoft record, Dewesoft's declared acquisition sample rate is
authoritative only after the selected-channel timestamp steps have been checked
against it. For the validated APE example this is exactly

[
f_s = 20000.0 mathrm{Hz}.
]

The downstream uniform time coordinate is therefore

[
t_n = t_0 + n/f_s.
]

The reader does not replace this metadata value with a frequency re-estimated
from floating-point timestamp differences.

## Calibration configuration

The same run YAML records four later calibration inputs:

[
C_{\mathrm{trans}},quad q,quad G,quad H_2.
]

They are run-specific configuration values, not universal APE constants.
Stage 1 does **not** apply them.

`second_harmonic_pd_response` is (H_2), the measured
photodiode/electronics response at (2f). For example, (H_2=0.7) means the
readout chain retains about 70% of the relevant (2f) response, so later
calibration must account for that attenuation explicitly.

The exact RMS/peak convention and equations producing (I_0(t)),
(eta_0(t)), (ho_X(t)), and (ho_Y(t)) belong to the separately
approved calibration stage, not to the raw reader.


## Streaming production reader

Production-size Dewesoft files are validated in bounded chunks rather than
loaded into one giant in-memory object. The default wrapper chunk size is
1,000,000 samples per channel.

For each chunk the reader verifies identical sample counts, synchronized
timestamps across all five channels, and agreement with Dewesoft's declared
sample rate. It also verifies continuity between the last timestamp of one
chunk and the first timestamp of the next.

Chunk boundaries are only a memory-management device. They do not remove,
insert, resample, reorder, interpolate, or otherwise alter samples.

`inspect_dxd` validates the complete file and returns metadata only.
`iter_ape_dxd_chunks` yields validated raw chunks for later pipeline stages.
The compatibility `read_ape_dxd` materializes all chunks and should not be
used for production-size acquisitions.
