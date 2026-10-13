# Ingestion-to-Spectrum Architecture

## Canonical APE stage 1

The new first pipeline boundary is intentionally independent of cleaning,
calibration, and spectral estimation:

```text
configured Dewesoft .dxd file
            |
            v
five explicit recorded channel names
X1f  Y1f  X2f  Y2f  cavity transmission
            |
            v
       read_ape_dxd
            |
            v
        APEDXDData
```

The reader preserves raw recorded values and timing metadata. It refuses
missing channels, duplicate mappings, unequal lengths, or unsynchronized
timestamps. It never trims, resamples, interpolates, filters, fills, cleans, or
calibrates data.

For uniformly sampled synchronized channels, the object stores (t_0),
(f_s), and (N) rather than a redundant full time array. If timestamps are
irregular, the actual relative timestamps are retained and no fake uniform
sample rate is assigned.

For production selection of the longest usable segment, the pipeline uses a
cavity-first ranked search.

```text
full declared-continuous acquisition
        |
        v
read cavity transmission only
        |
        v
0.05 s cavity mean/minimum
        |
        v
centered 60 s cavity q90
90% lower cut, 3 refinements
        |
        v
rank cavity-good segments by duration
        |
        v
read X2f, Y2f, and cavity only for the longest candidate
        |
        v
sample-level eta0
        |
        v
0.05 s eta0 mean/minimum
        |
        v
past-only accepted-history q90
80% eta0-minimum lower cut
        |
        v
longest surviving eta0-good subsegment
        |
        v
evaluate next cavity candidate only if it can still win
```

Eta0 protected-reference state is reset and bootstrapped at each cavity-good
segment start. It is not reset merely because a candidate crosses a declared-
continuous DXD file boundary. Failed eta0 blocks do not update the reference;
if the accepted-history window empties during a sustained eta0 failure, the
last healthy reference is held.

The cavity-only pass reads only the cavity transmission channel. Full
X1f/Y1f/X2f/Y2f/cavity data are read later only over candidate raw sample
ranges. This reduces I/O and avoids sample-level eta0 calibration across data
that cannot belong to the final longest segment.

The stopping rule is exact under the approved candidate-reset convention:
once the current best final-good duration is at least as long as the next
unprocessed cavity-good candidate, no later candidate can win because eta0
cleaning can only shorten it.

Raw science channels remain bounded-memory and are never interpolated,
resampled, smoothed, or gap-filled by this selection stage.

## Frozen calibrated HDF5 product

After the frozen segment is selected, production calibration reads only the
manifest-listed raw slices, in order, and applies
`ape_rms_quadrature_calibration_v1` chunk by chunk. The canonical calibrated
product is one HDF5 file with four aligned float64 datasets:

- `i0_w`
- `eta0`
- `rho_x`
- `rho_y`

The HDF5 root metadata records sample rate, sample count, frozen-manifest
SHA-256, source parts, calibration constants, calibration method, and code
commit. X and Y remain separate. This stage applies no filtering, detrending,
resampling, interpolation, smoothing, or X/Y combination.

Production writing is transactional: data are first written to a temporary
partial file; the final HDF5 path appears only after every manifest sample has
been calibrated, every calibrated sample is finite, the total length matches
the frozen manifest, and the product-level checks pass. Existing output files
are never overwritten.

## Legacy/reference ingestion and segment selection

The earlier chronology/quality-review architecture remains available for
cross-checking:

```text
explicit ordered source manifest / manual slices
                    |
                    v
reference Modules 2A + 2B + 3 (unchanged calls)
                    |
                    v
per-file quality masks and review artifacts
                    |
        split at invalid/gap/unverified boundary
                    |
                    v
longest globally continuous valid segment
                    |
                    v
immutable frozen validated-segment manifest
             /                         \
            v                           v
separate X voltage                 separate Y voltage
   -> FFT extension                   -> FFT extension
   -> C-LPSD backend                  -> C-LPSD backend
```

That path is not silently invoked by the canonical stage-1 reader.

Manual and automatic selection are interchangeable only at the frozen-manifest
interface. Both retain one channel per product. Combined plots, if approved,
are display-only references to separate method-specific artifacts.

For each diagnostic channel, legacy review plots are produced independently at
raw, cleaned/mask, and selected-segment stages. Rendering envelopes are not
analysis arrays; full-resolution focus plots are bounded by YAML
`zoom_max_samples`.
