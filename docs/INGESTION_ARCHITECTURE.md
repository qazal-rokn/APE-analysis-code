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

The approved downstream order is now:

```text
validated raw DXD chunks
        |
        v
0.05 s reduced records
        |
        v
cavity average -> centered 60 s q90 reference
cavity minimum -> 90% lower decision, 3 refinements
        |
        v
production cavity-good mask
        |
        v
sample-level eta0 from X2f/Y2f/cavity transmission
        |
        v
0.05 s eta0 mean + eta0 minimum
        |
        v
past-only 60 s q90 from previously accepted eta0 means
eta0 minimum -> 80% lower decision
failed eta0 blocks do not update reference
hold last healthy reference if accepted-history window empties
        |
        v
final good reduced-record mask
        |
        v
expand mask to original raw indices
        |
        v
full rho_X / rho_Y calibration only on final-good samples
```

The cavity quality reference is constructed across the full record at 0.05 s
resolution so its centered window does not reset at raw reader chunk
boundaries. The eta0 quality state is also continuous across declared
continuous file/chunk boundaries: its previous-60-s accepted history is not
reset there. The raw science channels remain bounded-memory and chunked.
Reduced statistics define quality only; they are not substituted for the raw
science samples used later in calibration and spectral estimation.

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
