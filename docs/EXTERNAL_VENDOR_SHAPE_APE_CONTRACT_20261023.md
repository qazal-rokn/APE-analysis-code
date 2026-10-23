# External scalar-field peak shape reused for APE (draft, 2026-10-23)

## Approval
Researcher explicitly requested unchanged scalar-vendor normalized signal shape and likelihood, with scalar LIGO detector response replaced by supplied APE QWP response. No new branch. Do not fabricate an intrinsic shape or resample without a separate approved response contract.

## Source behavior verified
The original external `max_lkl_executable.py` loads `peak_shape=np.load(peak_shape_file)` or defaults to `np.ones(1)`. Its `stats.log_likelihood` adds `peak_norm * peak_shape * mu` to the linear fitted background before calculating skew-normal log residual likelihood. The original detector-specific normalization is `stats.get_beta(freq_Hz, f_A_star)`; this must not be reused for APE polarization.

## APE implementation
`analysis/vendor_shape_ape.py` accepts the unmodified external .npy shape with sum=1 and the exact native C frequency-bin vector. It supplies `peak_norm[q,j] = (g_ref R_APE(f_j))^2 alpha_q/2`, where `R_APE` is the existing QWP nominal peak-amplitude response, and `mu=(g/g_ref)^2`. This preserves the external `background + mu*peak_norm*peak_shape` mathematical structure while replacing only detector response. Log observed/background PS are in rad².

## Crucial applicability condition
Length matching and sum-to-one alone do NOT show the external input shape includes native APE Constant-Q window leakage. The direct-shape adapter is physically appropriate only if the shape is defined as probability mass over the APE *observed estimator bins*, including window response. Otherwise compute the convolved physical template with `prepare_ape_likelihood` and the verified C kernel instead. Do not double-convolve and do not silently reinterpret intrinsic halo probability as measured-bin power. Sampling duration, coherence and X/Y response may also differ by detector.

## Evidence and pending gates
Focused contract tests were committed; no executable or genuine external peak-shape file was available for this run. Obtain the exact external `peak-shape-file` and its frequency definition and locate the archived PS-enabled C executable with independent binary SHA-256. Verify its native kernel, QWP RMS convention and X/Y signal fractions by injections. Original vendor/C, `main` and data unchanged.
