"""Bounded paper-2026 halo injection -> original native C -> original vendor.

Science:
 * Eq.(5) isotropic virialized Galactic 3D Gaussian, truncated at escape speed.
 * Fixed benchmark lab speed (Earth orbit not modeled) and nonrelativistic
   f_obs = f_rest * (1 + |v_gal-v_lab|**2/(2*c**2)).
 * One COMMON random-phase superposition of velocity modes produces the real
   X/Y time series. No independent per-bin signal powers.
 * APE's existing provisional QWP amplitude response, and an explicitly
   SURROGATE 90-degree X/Y phase relation; no complete physical phase model.
 * Explicit independent white readout noise. Real-data inference forbidden.

One fixed seed, one injected pair and its PAIRED noise-only control. The two
X/Y pairs are passed to the compiled native C backend and unchanged original
vendor, then the output is scored with the original vendor log likelihood
using a continuous halo's C-bin Kaiser response. No new C implementation and
no empirical K correction. This is an injection recovery test, NOT coverage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c

from physics.qwp_signal import qwp_amplitude_per_coupling, qwp_signal_power_template
from ape_c_lpsd.analysis.vendor_statistics import fit_window, factor_power_template, log_likelihood
from tools.run_ape_pipeline import execute

PAPER_DOI = "10.1038/s41598-025-33428-2"
F_REST = 5000.0
THETA_0 = 229000.0
SIGMA_3D_COMPONENT = THETA_0 * np.sqrt(1.5)
V_ESC = 528000.0
# Fixed solar/lab benchmark only, NOT an actual observation epoch.
V_LAB = 229000.0
SEED = 20261028
DURATION = 201.0
FS = 20000.0
N_BINS = 1301
PSLL = 238.13
Q = 1000000.0


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def sample_modes(rng, count, *, f_rest=F_REST):
    """Rejection-sample Eq.(5) speed distribution via 3D Galactic Gaussians."""
    if not (isinstance(count, int) and 64 <= count <= 8192 and f_rest > 0):
        raise ValueError("Choose 64..8192 modes and a positive rest frequency")
    velocities = []
    remaining = count
    while remaining:
        trial = rng.normal(0, SIGMA_3D_COMPONENT, (max(remaining*2, 128), 3))
        valid = trial[np.sum(trial**2, axis=1) <= V_ESC**2]
        accepted = valid[:remaining]
        velocities.append(accepted)
        remaining -= len(accepted)
    gal = np.concatenate(velocities)
    obs = gal - np.array([V_LAB, 0.0, 0.0])
    frequencies = f_rest*(1.0 + np.sum(obs**2, axis=1)/(2*c**2))
    phases = rng.uniform(0, 2*np.pi, count)
    indices = np.argsort(frequencies)
    return frequencies[indices], phases[indices]


def baseband_envelope(frequencies, phases, response_per_coupling, coupling, *,
                      f_rest, duration, envelope_rate=8., chunk_modes=128):
    """Continuous off-grid plane waves evaluated only on a cheap baseband grid.

    Linear interpolation of a <=0.02-Hz envelope sampled at 8 Hz introduces a
    negligible interpolation approximation (NOT an axion linewidth filter).
    """
    n = int(round(duration*envelope_rate)) + 1
    time = np.arange(n)/envelope_rate
    output = np.zeros(n, dtype=np.complex128)
    count = len(frequencies)
    for lo in range(0, count, chunk_modes):
        hi = min(lo+chunk_modes, count)
        angle = 2*np.pi*(frequencies[lo:hi, None]-f_rest)*time[None, :]
        angle += phases[lo:hi, None]
        output += np.sum(
            response_per_coupling[lo:hi, None]*np.exp(1j*angle), axis=0)
    output *= coupling/np.sqrt(count)
    return time, output


def produce_pair(out, *, config, count, seed, coupling, noise_rms, block_size=262144):
    """Stream paired injected/noise-only HDF5; no large RAM allocation."""
    if not (coupling > 0 and noise_rms > 0):
        raise ValueError("Positive benchmark injected coupling and noise RMS required")
    if not (0 < FS/2 and DURATION > 200.14):
        raise ValueError("Insufficient coherence-matched duration")
    rng = np.random.default_rng(seed)
    frequencies, phases = sample_modes(rng, count)
    phys = config["physics"]
    fractions = np.asarray(phys["quadrature_power_fraction"], float)
    if (fractions.shape != (2,) or not np.isfinite(fractions).all()
            or np.any(fractions < 0) or not np.isclose(fractions.sum(), 1.)):
        raise ValueError("Require explicit X/Y quadrature fractions summing to 1")
    response = qwp_amplitude_per_coupling(
        frequencies, finesse_s=phys["finesse_s"], tau_s=phys["tau_s"],
        tau1_s=phys["tau1_s"], rho_dm_gev4=phys["rho_dm_gev4"],
        unit_conversion=phys["unit_conversion"])
    slow_t, slow_field = baseband_envelope(
        frequencies, phases, response, coupling,
        f_rest=F_REST, duration=DURATION)
    n = round(FS*DURATION)
    paths = {key: out/(key+".h5") for key in ("injected", "noise_only")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite synthetic data")
    handles = {}
    try:
        for key, path in paths.items():
            f = h5py.File(path, "x")
            handles[key] = f
            for channel in ("rho_x", "rho_y"):
                f.create_dataset(channel, shape=(n,), dtype="f8",
                                 chunks=(min(block_size,n),))
            f.create_dataset("signal_frequency_hz", data=frequencies)
            f.create_dataset("vendor_shape_weights", data=np.full(count, 1/count))
            for attr, val in {
                "synthetic_only": True, "unit": "rad",
                "model": "stochastic_gaussian_vendor_shape",
                "source_paper_doi": PAPER_DOI, "seed": seed,
                "fs_hz": FS, "duration_s": DURATION,
                "rest_frequency_hz": F_REST,
                "injected_coupling_gev_inv": coupling if key == "injected" else 0.0,
                "qwp_phase_response": "envelope_only_90_degree_XY_surrogate",
                "halo_velocity_model": "Galactic_isotropic_3D_Gaussian_truncated_boosted_lab_fixed_no_orbit",
                "halo_theta0_m_s": THETA_0, "halo_vesc_m_s": V_ESC,
                "halo_lab_speed_m_s": V_LAB,
            }.items():
                f.attrs[attr] = val
        for start in range(0, n, block_size):
            stop = min(start+block_size, n)
            t = np.arange(start,stop,dtype=float)/FS
            xr = np.interp(t,slow_t,slow_field.real)
            yi = np.interp(t,slow_t,slow_field.imag)
            carrier = 2*np.pi*F_REST*t
            cos = np.cos(carrier)
            sin = np.sin(carrier)
            # Real and imaginary parts of same analytic field, 90-degree
            # X/Y surrogate. Real signal powers match QWP T=0.5*g^2*alpha*R^2.
            sx = np.sqrt(fractions[0])*(xr*cos-yi*sin)
            sy = np.sqrt(fractions[1])*(xr*sin+yi*cos)
            nx = rng.normal(0.,noise_rms,len(t))
            ny = rng.normal(0.,noise_rms,len(t))
            handles["noise_only"]["rho_x"][start:stop] = nx
            handles["noise_only"]["rho_y"][start:stop] = ny
            handles["injected"]["rho_x"][start:stop] = nx+sx
            handles["injected"]["rho_y"][start:stop] = ny+sy
    finally:
        for h in handles.values():
            h.close()
    return paths, frequencies


def window_power_kernel(output_frequency_hz, mode_frequency_hz, *, q, psll):
    """Native C Kaiser normalized coherent-tone power response, no sum scaling."""
    f = np.asarray(output_frequency_hz, float)
    modes = np.asarray(mode_frequency_hz, float)
    if not np.isfinite(f).all() or not np.isfinite(modes).all() or np.any(f <= 0):
        raise ValueError("Invalid frequencies")
    x = psll/100.
    alpha = ((.0889732*x-.493285)*x+4.71469)*x-.0821377
    u = (modes[None,:]-f[:,None])*q/f[:,None]
    z = alpha*alpha-u*u
    numer = np.empty_like(u)
    inside = z >= 0
    arg = np.pi*np.sqrt(np.maximum(z,0))
    numer[inside] = np.where(arg[inside] == 0, 1., np.sinh(arg[inside])/np.where(arg[inside] == 0, 1., arg[inside]))
    arg = np.pi*np.sqrt(np.maximum(-z,0))
    numer[~inside] = np.sin(arg[~inside])/np.where(arg[~inside]==0,1,arg[~inside])
    denom = np.sinh(np.pi*alpha)/(np.pi*alpha)
    return (numer/denom)**2


def analyze_native(out, config, frequencies, injected_coupling):
    """Original vendor background + log-likelihood, native C PS unchanged."""
    vroot = config["external_vendor_root"]
    spectra = {}
    fits = {}
    for key in ("noise_only","injected"):
        directory = out/("ape-pipeline-"+key)
        pair = []
        for channel in ("rho_x","rho_y"):
            path=directory/("ape-ps-"+channel)/(channel+".joined.ps")
            data=np.loadtxt(path,ndmin=2)
            if data.shape != (N_BINS,2):
                raise ValueError("Unexpected native C shape")
            pair.append(data)
        if not np.array_equal(pair[0][:,0], pair[1][:,0]):
            raise ValueError("Native C X/Y frequencies differ")
        if np.any(np.stack([p[:,1] for p in pair]) <= 0):
            raise ValueError("Nonpositive native C power")
        spectra[key] = np.stack([p[:,1] for p in pair])
        fits[key] = [fit_window(pair[0][:,0], p[:,1], vroot) for p in pair]
        if any(not x["vendor_quality_pass"] for x in fits[key]):
            raise RuntimeError("Original vendor background fit failed quality criterion")
    native_f = pair[0][:,0]
    ratios = native_f[1:]/native_f[:-1]
    q_actual = float(1/np.median(ratios-1))
    if abs(q_actual/Q-1)>1e-3:
        raise ValueError("Actual native C grid is NOT paper Q=1e6")
    W = window_power_kernel(native_f, frequencies, q=q_actual, psll=PSLL)
    alpha = np.asarray(config["physics"]["quadrature_power_fraction"],float)
    physics = config["physics"]
    ref = float(physics["reference_coupling_gev_inv"])
    template = qwp_signal_power_template(
        signal_frequency_hz=frequencies,
        normalized_halo_weights=np.full(len(frequencies),1/len(frequencies)),
        window_power_response=np.broadcast_to(W[None,:,:],(2,*W.shape)),
        quadrature_power_fraction=alpha, reference_coupling_gev_inv=ref,
        finesse_s=physics["finesse_s"], tau_s=physics["tau_s"],
        tau1_s=physics["tau1_s"], rho_dm_gev4=physics["rho_dm_gev4"],
        unit_conversion=physics["unit_conversion"])
    projected = np.sum(template,axis=0)
    if not np.isfinite(projected).all() or projected.max()<=0:
        raise ValueError("Invalid signal PS template")
    active = np.flatnonzero(projected >= projected.max()*1e-4)
    # This selection is derived from actual modeled signal, NOT a fixed 7 bins.
    lo,hi = int(active[0]),int(active[-1]+1)
    if hi-lo>80 or lo<0 or hi>N_BINS:
        raise RuntimeError("Modeled native signal too broad or invalid")
    slice_ = slice(lo,hi)
    norm, shape = factor_power_template(template[:,slice_])
    coupling_grid = np.array([0., .25, .5, .75, 1., 1.25, 1.5, 2.])*injected_coupling
    report = {}
    for key in ("noise_only","injected"):
        P = spectra[key][:,slice_]
        B = np.stack([fit["background_psd"][slice_] for fit in fits[key]])
        pars = np.stack([fit["skewnormal_loc_scale_alpha"] for fit in fits[key]])
        ll = np.array([log_likelihood(
            (float(g)/ref)**2, np.log(P),np.log(B),
            norm,shape,pars,vroot) for g in coupling_grid])
        if not np.isfinite(ll).all():
            raise RuntimeError("Original vendor returned nonfinite likelihood")
        best = int(np.argmax(ll))
        residual = np.sum(np.maximum(P-B,0),axis=0)
        peak_idx = lo+int(np.argmax(residual))
        report[key] = {
            "vendor_background_quality": [float(f["vendor_reduced_chi2"]) for f in fits[key]],
            "likelihood_grid": ll.tolist(),
            "relative_to_null_log_likelihood": (ll-ll[0]).tolist(),
            "best_grid_coupling_gev_inv": float(coupling_grid[best]),
            "max_excess_power_frequency_hz": float(native_f[peak_idx]),
            "local_composite_diagnostic_q0_not_significance": (
                float(max(0.,2*(ll[best]-ll[0])))),
        }
    return {
        "native_C_grid_Q": q_actual,
        "native_C_selected_frequency_bins_hz": native_f[slice_].tolist(),
        "signal_selected_bins_count_not_fixed_seven": int(hi-lo),
        "template_per_XY_PS_rad2_at_g_ref": template[:,slice_].tolist(),
        "coupling_grid_gev_inv": coupling_grid.tolist(),
        "halo_min_frequency_hz": float(frequencies.min()),
        "halo_max_frequency_hz": float(frequencies.max()),
        "halo_frequency_05_95_quantiles_hz": np.quantile(
            frequencies,[.05,.95]).tolist(),
        "realization_results": report,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--existing-config",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--coupling-gev-inv",type=float,required=True)
    p.add_argument("--noise-rms-rad",type=float,required=True)
    p.add_argument("--modes",type=int,default=1024)
    p.add_argument("--seed",type=int,default=SEED)
    p.add_argument("--run-c",action="store_true")
    args=p.parse_args()
    existing=Path(args.existing_config).expanduser().resolve(strict=True)
    config=json.loads(existing.read_text())
    dest=Path(args.output_dir).expanduser().resolve()
    if dest.exists():
        raise FileExistsError("Output dir must be new; no overwrite")
    if not dest.name.startswith("ape-paper-axion-"):
        raise ValueError("Require new ape-paper-axion-* run directory")
    for key in ("ps_executable","ps_sha256","external_vendor_root","physics","c_settings"):
        if key not in config:
            raise ValueError("Missing established approved configuration: "+key)
    exe=Path(config["ps_executable"]).expanduser().resolve(strict=True)
    if sha256(exe).lower()!=str(config["ps_sha256"]).lower():
        raise ValueError("Original C binary provenance mismatch")
    native_setting=dict(config["c_settings"])
    required = {"fs_hz","duration_s","f_min_hz","f_max_hz","grid_bins",
                "overlap_percent","window","psll_db","regression",
                "input_scale","epsilon","memory_power","kernel_threshold",
                "a_option","batch_count"}
    if set(native_setting)!=required:
        raise ValueError("C settings incomplete")
    if not (int(native_setting["window"]) == -2 and
            np.isclose(native_setting["psll_db"],PSLL) and
            int(native_setting["overlap_percent"]) == 20 and
            int(native_setting["a_option"]) == 2):
        raise ValueError("Require native C method-2 Kaiser PSLL 238.13 and 20 percent overlap")
    if not np.isfinite([args.coupling_gev_inv,args.noise_rms_rad]).all():
        raise ValueError("Nonfinite amplitude or noise")
    # Centered bin j=650 is EXACT f_rest; ratio yields exact Q=1e6.
    step=np.log1p(1/Q)
    native_setting.update(
        fs_hz=FS, duration_s=DURATION, grid_bins=N_BINS,
        f_min_hz=F_REST*np.exp(-650*step),
        f_max_hz=F_REST*np.exp(650*step))
    # Run at most one C injection and one paired null control.
    dest.mkdir(parents=True,exist_ok=False)
    paths,modes = produce_pair(
        dest,config=config,count=args.modes,seed=args.seed,
        coupling=args.coupling_gev_inv,noise_rms=args.noise_rms_rad)
    report = {
        "status": "paper_2026_halo_axion_like_time_domain_injection_created",
        "paper_doi": PAPER_DOI,
        "science_assumptions": {
            "Galactic_isotropic_halo_dispersion_component_m_s": float(SIGMA_3D_COMPONENT),
            "theta0_m_s": THETA_0, "escape_speed_m_s": V_ESC,
            "fixed_lab_speed_m_s": V_LAB, "earth_orbit_included": False,
            "Q_C": Q, "rest_frequency_hz": F_REST,
            "coherence_scale_approx_s": 1/(F_REST*1e-6),
            "halo_mode_count": args.modes, "field_seed": args.seed,
            "X_Y_mapping": "90_degree_surrogate_same_field_not_validated_APE_optical_phase",
            "noise": "two_independent_white_Gaussian_X_Y_paired_between_injected_and_null",
            "signal_stochasticity": "single_continuous_superposition_of_random_phase_velocity_modes",
        },
        "injection": {
            "g_injected_gev_inv": args.coupling_gev_inv,
            "white_noise_rms_rad": args.noise_rms_rad,
            "duration_s": DURATION,"sample_rate_hz":FS,
            "frequency_min_max_hz": [float(modes.min()),float(modes.max())],
            "frequency_quantile_05_95_hz": np.quantile(modes,[.05,.95]).tolist(),
        },
        "provenance": {
            "input_config_sha256": sha256(existing),
            "original_C_SHA256": sha256(exe),
            "injected_hdf5_sha256": sha256(paths["injected"]),
            "null_hdf5_sha256": sha256(paths["noise_only"]),
        },
        "native_c_executed": False,
        "original_vendor_executed": False,
        "physical_APE_phase_and_coverage_validated": False,
    }
    if args.run_c:
        for key in ("noise_only","injected"):
            cfg=dict(config)
            cfg.update(
                synthetic_only=True,execute_c=True,
                input_kind="synthetic_halo_file",require_stochastic_halo=True,
                input_hdf5=str(paths[key]),run_dir=str(dest/("ape-pipeline-"+key)),
                c_settings=native_setting,peak_shape_file=None,
                statistical_mode="finite_grid",candidate_start_indices=[
                    642,646,650,654,658,662],
                coupling_grid_gev_inv=[
                    0., .25*args.coupling_gev_inv,
                    .5*args.coupling_gev_inv,
                    args.coupling_gev_inv,
                    1.5*args.coupling_gev_inv],
                n_effective_trials=None,local_sigma_threshold=3.0)
            # All C and original vendor processing is delegated unchanged.
            completed=execute(cfg)
            report["provenance"][key+"_native_manifest"]=str(
                dest/("ape-pipeline-"+key)/"manifest.json")
            report["provenance"][key+"_candidate_count"]=len(
                completed["candidates"]["candidates"])
        report["native_likelihood_comparison"]=analyze_native(
            dest,config,modes,args.coupling_gev_inv)
        report["native_c_executed"]=True
        report["original_vendor_executed"]=True
        report["status"]="PAPER_2026_AXION_LIKE_NATIVE_C_AND_VENDOR_RECOVERY_DIAGNOSTIC_COMPLETE"
    with (dest/"axion_like_recovery_report.json").open("x") as handle:
        json.dump(report,handle,indent=2)
    print(report["status"])
    print("NATIVE_C_EXECUTED",report["native_c_executed"])
    print("ORIGINAL_VENDOR_EXECUTED",report["original_vendor_executed"])
    print("HALO_FREQUENCY_RANGE_HZ",report["injection"]["frequency_min_max_hz"])
    if args.run_c:
        r=report["native_likelihood_comparison"]
        print("NATIVE_C_Q",r["native_C_grid_Q"])
        print("PHYSICALLY_SELECTED_BIN_COUNT",r["signal_selected_bins_count_not_fixed_seven"])
        for key in ("noise_only","injected"):
            x=r["realization_results"][key]
            print(key,"best_grid_g",x["best_grid_coupling_gev_inv"],
                  "frequency_max_excess",x["max_excess_power_frequency_hz"])
    print("REPORT",dest/"axion_like_recovery_report.json")
    print("NO_NEW_NORMALIZATION_NO_PHYSICAL_LIMIT_OR_SIGNIFICANCE")


if __name__ == "__main__":
    main()
