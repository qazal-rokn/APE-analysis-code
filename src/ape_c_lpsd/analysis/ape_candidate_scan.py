"""Provisional APE candidate likelihood scan using unchanged external vendor stats.

Accepts already calibrated, fitted candidate windows. No data simulation or
real-data authorization; finite grid search is a diagnostic, not coverage.
"""
import numpy as np
from ape_c_lpsd.analysis.vendor_shape_ape import prepare_vendor_shape_ape_likelihood
from ape_c_lpsd.analysis.ape_candidates import select_candidates
from ape_c_lpsd.analysis.vendor_statistics import log_likelihood


def scan_candidate_windows(*, candidate_frequencies_hz, observed_ps,
                           background_ps, skewnormal_parameters,
                           frequency_windows_hz, vendor_peak_shape,
                           quadrature_power_fraction, reference_coupling_gev_inv,
                           coupling_grid_gev_inv, external_root,
                           local_sigma_threshold, finesse_s, tau_s,
                           tau1_s, rho_dm_gev4, unit_conversion,
                           n_effective_trials=None):
    """Scan user-selected frequency windows; first axis is candidate index.

    observed/background shape (C,M,N); frequency_windows shape (C,N);
    skewnormal shape (C,M,3). No interpolation or independent noise claim.
    """
    centers=np.asarray(candidate_frequencies_hz,dtype=float)
    frequency=np.asarray(frequency_windows_hz,dtype=float)
    obs=np.asarray(observed_ps,dtype=float)
    bg=np.asarray(background_ps,dtype=float)
    pars=np.asarray(skewnormal_parameters,dtype=float)
    grid=np.asarray(coupling_grid_gev_inv,dtype=float)
    if (centers.ndim!=1 or not centers.size or not np.isfinite(centers).all()
            or np.any(np.diff(centers)<=0) or frequency.ndim!=2
            or frequency.shape[0]!=len(centers)
            or obs.ndim!=3 or bg.shape!=obs.shape
            or obs.shape[0]!=len(centers) or obs.shape[2]!=frequency.shape[1]
            or pars.shape!=(obs.shape[0],obs.shape[1],3)
            or grid.ndim!=1 or grid.size<2 or not np.isfinite(grid).all()
            or np.any(grid<0) or np.any(np.diff(grid)<=0) or grid[0]!=0):
        raise ValueError("Require aligned (C,M,N) inputs and sorted coupling grid starting at zero")
    values=[]
    for i in range(len(centers)):
        prepared=prepare_vendor_shape_ape_likelihood(
            observed_ps=obs[i],background_ps=bg[i],
            skewnormal_parameters=pars[i],frequency_bins_hz=frequency[i],
            vendor_peak_shape=vendor_peak_shape,
            quadrature_power_fraction=quadrature_power_fraction,
            reference_coupling_gev_inv=reference_coupling_gev_inv,
            finesse_s=finesse_s,tau_s=tau_s,tau1_s=tau1_s,
            rho_dm_gev4=rho_dm_gev4,unit_conversion=unit_conversion)
        val=np.array([log_likelihood(
            (float(g)/reference_coupling_gev_inv)**2,prepared["Y"],
            prepared["bkg"],prepared["peak_norm"],prepared["peak_shape"],
            prepared["model_args"],external_root) for g in grid])
        if not np.isfinite(val).all():
            raise ValueError("Vendor returned non-finite likelihood scan")
        values.append(val)
    values=np.asarray(values)
    best=np.argmax(values,axis=1)
    return select_candidates(centers,values[:,0],values[np.arange(len(best)),best],
                             (grid[best]/reference_coupling_gev_inv)**2,
                             local_sigma_threshold=local_sigma_threshold,
                             n_effective_trials=n_effective_trials) | dict(
        scan_method="finite_grid_provisional_vendor_likelihood",
        coupling_grid_gev_inv=grid.tolist(),
        best_fit_coupling_gev_inv=grid[best].tolist(),
        log_likelihood_null=values[:,0].tolist(),
        log_likelihood_best=values[np.arange(len(best)),best].tolist())
