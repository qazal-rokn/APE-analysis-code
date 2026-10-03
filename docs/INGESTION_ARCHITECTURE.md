# Ingestion-to-Spectrum Architecture

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

Manual and automatic selection are interchangeable only at the frozen-manifest
interface. Both retain one channel per product. Combined plots, if approved,
are display-only references to separate method-specific artifacts.

For each diagnostic channel, review plots are produced independently at raw,
cleaned/mask, and selected-segment stages. Rendering envelopes are not analysis
arrays; full-resolution focus plots are bounded by YAML `zoom_max_samples`.
