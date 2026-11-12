# Canonical APE Raw-Data Contract

The first APE pipeline stage is deliberately narrow:

```text
Dewesoft .dxd
    |
    v
five explicitly named recorded channels
    |
    v
bounded synchronized chunks
```

No calibration, filtering, cleaning, interpolation, resampling, trimming, gap
filling, segment selection, spectral estimation, or X/Y combination is allowed
inside this stage.

## Required channels

The run configuration must map exactly five APE roles to five distinct Dewesoft
channel names:

- `x_1f`: recorded first-harmonic X quadrature;
- `y_1f`: recorded first-harmonic Y quadrature;
- `x_2f`: recorded second-harmonic X quadrature;
- `y_2f`: recorded second-harmonic Y quadrature;
- `cavity_transmission`: recorded cavity-transmission monitor.

The reader never guesses these names.

## Bounded-memory production interface

Large APE files are not materialized as one giant in-memory object.

`iter_ape_dxd_chunks(...)` reads the five synchronized raw channels in bounded
blocks. The default block contains 1,000,000 samples per channel. Each returned
`APEDXDChunk` contains:

- source file;
- raw start and stop sample indices;
- exact Dewesoft sample rate;
- the five raw channel arrays for that interval;
- the first and last Dewesoft timestamps for that interval.

The reader validates every chunk before yielding it.

`inspect_dxd(...)` exhausts the same chunked path across the entire file and
returns `APEDXDInspection`, which contains metadata only. The full recording
is therefore never retained in RAM by the inspection wrapper.

The small-input `APEDXDData` object remains only for controlled synthetic
helpers and tests. It is not the production path for multi-hour recordings.

## Validation rules

Before and during streaming the reader verifies:

1. all five configured channel names exist and are distinct;
2. all five channels are scalar channels;
3. all five channels contain the same total number of samples;
4. every chunk contains finite recorded values and finite timestamps;
5. timestamps agree across all five channels;
6. timestamp spacing agrees with Dewesoft's declared sample rate;
7. the first timestamp of each chunk follows continuously from the last
   timestamp of the previous chunk.

If any rule fails, ingestion stops. The reader does not truncate, align,
resample, interpolate, repair, or fill data.

## Timing rule

For a regular Dewesoft record, the downstream time coordinate is defined by
the file start and Dewesoft's declared sampling rate after the timestamp checks:

[
t_n = t_0 + rac{n}{f_s}.
]

For the validated real APE smoke-test file, Dewesoft declares

[
f_s = 20000.0 mathrm{Hz}.
]

The timestamp sequence is used as an independent consistency check, not as a
reason to replace the declared sample rate with a floating-point estimate such
as 19999.999959 Hz.

## Memory behavior

For a three-hour, 20 kHz recording, retaining five complete float64 channels
would already require about 8.64 GB before temporary timestamp buffers and
other overhead. The chunked interface avoids that scaling by keeping only one
bounded interval in memory at a time.

Chunking is a transport/memory strategy only. It does not change sample values,
sample order, timing, units, or channel definitions.

## Calibration configuration

The same run YAML records four later calibration inputs:

[
C_{mathrm{trans}},quad q,quad G,quad H_2 .
]

They are run-specific configuration values, not universal APE constants.
Stage 1 does **not** apply them.

`second_harmonic_pd_response` is (H_2), the measured
photodiode/electronics response at (2f). For example, (H_2=0.7) means the
readout chain retains about 70% of the relevant (2f) response, so the later
calibration must account for that attenuation explicitly.

The exact RMS/peak convention and the equations producing (I_0(t)),
(eta_0(t)), (ho_X(t)), and (ho_Y(t)) belong to the separately
approved calibration stage, not to the reader.
