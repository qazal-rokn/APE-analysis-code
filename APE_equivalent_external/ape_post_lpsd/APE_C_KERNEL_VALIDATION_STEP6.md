# Step 6: Verify the APE native-C PSD response kernel, without changing C

## Source and motivation

The unchanged upstream C source `LPSD/IO.c` maps format column `D` to
`data->psd` (one-sided PSD), while `S` denotes power spectrum. The
APE method-2 configuration uses `FORMAT "fDNbri"`, so its second output
column is **PSD in rad²/Hz**. The Step-5 APE halo/QWP model requires a
native **PSD response per unit time-domain signal variance**, rather than an
analytically guessed Kaiser response or an unrelated PS kernel.

No original C, vendor statistics, APE preprocessing or calibrated real data
are modified by this step.

## Implemented

One new APE-only module, `ape_post_lpsd/validate_native_c_kernel.py`:

1. `make_tone_fixture(...)` writes a non-overwriting synthetic HDF5
   `/rho_x` **or** `/rho_y` sinusoid of known *peak* radian amplitude.
   No randomized axion phase, stochastic dark-matter signal or real-data input.
2. Run the **existing unmodified C method-3 + method-2 workflow outside this
   module**, for each tone, separately for X and Y, at TWO peak amplitudes
   and identical phase/duration/C settings except input and output paths.
   The module deliberately does not run C or substitute a Python estimator.
3. `validate_c_kernel(...)` streams original C outputs, requiring the exact
   native frequency grid, matching realized `N_avg` and `N_kernel`, and
   consistent two-amplitude scaling. Accept only explicitly declared raw
   six-column `fDNbri` or previously validated combined four-column format.
4. The measured frequency response is defined with **no extra correction**:

       K_D[m,j,k] = D_low[m,j,k] / (A_peak,low²/2)  [Hz^-1]

   The second amplitude is **only a consistency check**, not averaged into
   the kernel or used for empirical rescaling. Failing bins are rejected,
   not altered. The 2% default tolerance is a *synthetic validation setting*;
   its adequacy for production must be justified independently.
5. `export_kernel(...)` saves the kernel as an NPZ with fields
   `native_frequency_hz`, `output_bin_indices`,
   `halo_frequency_hz`, `native_c_psd_kernel_per_variance_hz`,
   and `metadata_json`. These arrays can be supplied **unchanged** to
   Step 5 `build_response_psd(...)`. File hashes and scaling diagnostics
   are recorded, and overwrites are prohibited.

## Execution instructions for a controlled local tone run

Run from `APE_equivalent_external/` in an installed NumPy + h5py environment.

1. Generate **separate** X/Y tone HDF5 inputs at two known amplitudes,
   e.g. from Python:

       from ape_post_lpsd.validate_native_c_kernel import make_tone_fixture
       make_tone_fixture("/new/run/tone_X_low.h5",
           channel="X", frequency_hz=5000., peak_amplitude_rad=1e-7,
           sample_rate_hz=20000., duration_s=201.)
       make_tone_fixture("/new/run/tone_X_high.h5",
           channel="X", frequency_hz=5000., peak_amplitude_rad=2e-7,
           sample_rate_hz=20000., duration_s=201.)

   Repeat for Y and every desired tone frequency. The 201 s example is a
   *test proposal*, not a guarantee of sufficient C coherence/window support
   for every frequency or instrument profile.
2. Use the **unchanged C executable** with the existing APE planner's
   historically recorded method-3/method-2 options, separately for each
   synthetic tone and quadrature. Preserve the command line, config,
   compiled binary SHA-256, input HDF5 hash, and all output files.
   Verify matching C flags, FFT/window normalization, phase, and duration.
   The present Python validator **cannot independently prove** the commands
   actually used merely from a supplied declaration.
3. Create a small JSON campaign using ABSOLUTE paths to the **C outputs**:

       {
         "native_spectrum_hdf5": "/path/to/APE_T0_T1_X.h5",
         "selected_indices": [100, 101, 102],
         "mode_frequencies_hz": [5000.0],
         "peak_amplitudes_rad": [1e-7, 2e-7],
         "columns": 6,
         "c_outputs": {
           "X": [["/run/X_low.lpsd", "/run/X_high.lpsd"]],
           "Y": [["/run/Y_low.lpsd", "/run/Y_high.lpsd"]]
         },
         "c_configuration_source": "/run/native_c_commands_and_settings.txt",
         "c_executable_sha256": "REPLACE_WITH_ACTUAL_64_HEX_SHA256"
       }

   Indices must refer to **real native C output bins**, and the chosen
   frequencies must lie within the original C output grid. Replace all
   placeholders; do not copy the example's indices blindly. The helper reads
   the possibly multi-million-bin native frequency array from the Step-1 HDF5
   instead of embedding it in JSON.

4. Validate without altering the vendor estimator:

       python -m ape_post_lpsd.validate_native_c_kernel \
         --campaign-json /new/run/campaign.json \
         --output-npz /new/run/verified_tone_kernel.npz

   A successful run proves *input consistency and amplitude-squared
   scaling for those supplied output files*. It does **not** independently
   prove the provenance of the native C executable or calibration records.

## Focused tests

       python -m unittest discover -s ape_post_lpsd/tests \
         -p 'test_validate_native_c_kernel.py' -v

Ten synthetic tests passed locally, including both C output formats, input
grid/metadata rejection, quadratic response, deterministic HDF5 tone creation,
HDF5 campaign loading, and Step-5-compatible NPZ export. All tests used
synthetic text spectra, NOT a compiled-C scientific validation run.

## Scientific status and next gate

**STATUS: VALIDATOR IMPLEMENTED, NATIVE C KERNEL NOT YET MEASURED.**
Before authorizing physical APE coupling inference, run controlled synthetic
tones through the real compiled C estimator at approved settings; verify its
binary and command provenance; check the PSD kernel and window support;
and independently validate cavity/two-QWP amplitude, phase response,
laboratory halo and X/Y/bin covariance assumptions.

The original vendor C, background, likelihood, optimizer, candidate and
upper-limit code remain byte-for-byte unchanged.
