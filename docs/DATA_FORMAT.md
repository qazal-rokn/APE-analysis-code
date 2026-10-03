# Canonical APE Raw-Data Contract

The first APE pipeline stage is deliberately narrow:

```text
Dewesoft .dxd
    |
    v
five explicitly named recorded channels
    |
    v
APEDXDData
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

## Output object

`APEDXDData` contains:

- source path;
- the five raw one-dimensional arrays;
- selected-channel names and recorded units;
- sample count;
- start-store time reported by Dewesoft, without inventing a timezone;
- sample rate when the selected timestamps are uniformly spaced;
- explicit relative timestamps if spacing is irregular;
- Dewesoft metadata used to audit the read.

If selected channel lengths or timestamps disagree, ingestion fails. The reader
does not truncate, align, resample, or repair them.

## Timing rule

For a regular record, the downstream time coordinate is defined by the file
start and the measured selected-channel sampling rate:

[
t_n = t_0 + n/f_s .
]

There is therefore no reason to retain another full time array for an ordinary
uniform record. If timing is not uniform, the actual relative timestamps are
preserved and `sample_rate_hz` is left unset rather than pretending that the
record is uniform.

## Calibration configuration

The same run YAML records four later calibration inputs:

[
C_{\mathrm{trans}},\quad q,\quad G,\quad H_2 .
]

They are run-specific configuration values, not universal APE constants.
Stage 1 does **not** apply them.

`second_harmonic_pd_response` is (H_2), the measured photodiode/electronics
response at (2f). For example, (H_2=0.7) means that the readout chain
retains about 70% of the relevant (2f) response, so the later calibration
must account for that attenuation explicitly.

The exact RMS/peak convention and the equations producing (I_0(t)),
(eta_0(t)), (ho_X(t)), and (ho_Y(t)) belong to the separately
approved calibration stage, not to the reader.
