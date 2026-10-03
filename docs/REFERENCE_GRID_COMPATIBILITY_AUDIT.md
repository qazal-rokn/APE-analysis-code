# Reference Grid Compatibility and Feasibility Audit

Date: 2026-09-25

## Input identity

Audited input: `reference_grid/final_native_65Hz_9kHz_exact_ccore_plan.npz`.
SHA-256 verified: `17d249bf18b98be6df3813adb441598b1b7931287fa7c383a24e0a46b7ced99f`.
The file was loaded with `allow_pickle=False`.

## Stored plan validation

- Arrays: `frequency_hz`, `r`, `m`, `L`, and `chunk_bounds`.
- `frequency_hz`, `r`, and `m` are finite positive float64 arrays of 5,025,124 elements.
- `L` is a positive int64 array of the same length.
- `frequency_hz` is strictly increasing, from 65.00004858924416 Hz to 8999.992444832225 Hz.
- `chunk_bounds` is int64, starts at 0, ends at 5,025,124, and its 130 positive intervals form an exact non-overlapping partition.
- The stored identity `L = round(20000/r)` holds exactly.
- With N=308,717,000 and fs=20,000 Hz, maximum stored L is 307,692,078 at node 0 (65.00004858924416 Hz). No stored node has L>N.

## Unchanged C CLI compatibility

The inspected `lpsd-exec` CLI accepts scalar `-s fmin`, `-t fmax`, `-J Jdes`, `-n nspec`, and `-N iter`; it does not accept an external frequency array or per-node L/m arrays. `calc_params` in `lpsd.c` generates `f_j=fmin*exp(j*log(fmax/fmin)/(Jdes-1))`. Batches select consecutive indices of this generated grid and therefore cannot substitute an arbitrary external grid.

Using the NPZ endpoints and Jdes=5,025,124 in that exact C formula gives one matching endpoint node only (node 0). Nodes 1 through 5,025,123 mismatch; the first difference is 1.006945339554477e-06 Hz at node 1 and the maximum difference is 62.875138183886065 Hz at node 4,015,595.

For that different C-generated grid, the scalar relative step is 9.811886034771788e-07 and the C LPSD DFT requirement is 313,591,165 samples at node 0; 15,966 generated nodes require more than N. This is not a failure of the supplied plan: its stored L values all fit. It demonstrates that regenerating a superficially similar C grid is not an acceptable substitute.

The source exposes `calculateSpectrum(tCFG*, tDATA*)`, `calculate_lpsd`, and constant-Q routines, all driven by tCFG scalar grid parameters and internally generated `data->fspec`; no existing C-core public interface accepts the NPZ node/L arrays.

## Conclusion

The supplied reference grid is valid and every stored required window fits the approved segment. It is incompatible with the unchanged `lpsd-exec` executable for exact-grid calculation: neither one invocation, non-overlapping `-n`/`-N` batches, nor an existing C-core API can consume or reproduce the node array exactly. No real signal data were read, no HDF5 was written, and no C executable was invoked during this audit.
