#!/usr/bin/env python3
"""APE one-pass, FAST combined likelihood validation at saved 5/8 kHz cases.

NO change to original vendor code or real data. One invocation checks:
  (1) halo -> original QWP -> C-D PSD forward-model numerical identity,
  (2) real 4.28h X/Y native grids, PSD and original vendor background fits,
  (3) measured local bin/X-Y correlation,
  (4) paired, correlated-noise null and additive-signal likelihood separation,
      using the UNCHANGED vendor stats.log_likelihood,
  (5) a few original vendor calc_max_lkl optimizer injection recoveries.

Monte Carlo is an EMPIRICAL COPULA SURROGATE: contiguous observed high-pass
logPSD neighborhoods preserve sampled rank dependence; marginal ranks are
mapped to the fitted vendor skew-normal residual distribution.  This DOES
NOT reproduce the C time-series estimator, axion stochastic coherence,
full-band trial factors, genuine coverage or physical g upper limits.
X/Y each use a separate counterfactual alpha=1; NEVER combine them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import rankdata, skewnorm

from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model
from scalar_dark_matter_lpsd.lvk.stats import log_likelihood
from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl
from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def jf(x):
    z = float(x)
    return z if np.isfinite(z) else None


def acf(z: np.ndarray, maxlag: int = 20):
    z = np.asarray(z, float)
    z = z-z.mean()
    out = np.ones(maxlag+1)
    for k in range(1, maxlag+1):
        out[k] = float(np.corrcoef(z[:-k], z[k:])[0, 1])
    return out


def inflation_ratio(T, B, rho):
    w = T/B
    iid = float(np.dot(w, w))
    cov = iid
    for k in range(1, min(len(rho), len(w))):
        cov += 2.0*float(rho[k])*float(np.dot(w[:-k], w[k:]))
    if iid <= 0 or cov <= 0:
        raise ValueError('Nonpositive weighted variance')
    return float(np.sqrt(cov/iid))


def vendor_fit(background, spectrum, freq, chi_max):
    table = background.get_df(background.get_label(filepath=str(spectrum)))
    fields = {'fmin', 'fmax', 'x_knots', 'y_knots', 'alpha_skew', 'loc_skew', 'sigma_skew', 'chi_sqr'}
    if not fields.issubset(table.columns):
        raise ValueError('Fitted vendor background missing required columns')
    same = table[(table['fmin'] <= float(freq[0])) & (table['fmax'] >= float(freq[-1]))]
    if same.empty:
        raise ValueError('81-bin band crosses vendor background segment boundary')
    fit = same.iloc[0]
    chi = float(fit['chi_sqr'])
    if not (np.isfinite(chi) and chi < chi_max):
        raise ValueError(f'Background chi-square {chi} failed < {chi_max}')
    knots = np.r_[np.asarray(fit['x_knots'], float), np.asarray(fit['y_knots'], float)]
    logB = np.asarray(bkg_model(np.log(freq), knots), float).reshape(-1)
    pars = np.array([[fit['loc_skew'], fit['sigma_skew'], fit['alpha_skew']]], float)
    if (logB.shape != freq.shape or not np.isfinite(logB).all() or
        not np.isfinite(pars).all() or pars[0,1] <= 0):
        raise ValueError('Invalid vendor background spline/skew-normal fit')
    B = np.exp(logB)
    if not np.isfinite(B).all() or np.any(B <= 0):
        raise ValueError('Invalid background PSD')
    return logB, B, pars, chi


def fit_original(Y, logB, T, pars, seeds):
    """Invoke UNMODIFIED signed-mu vendor optimizer for one 81-bin sample."""
    Y = np.asarray(Y, float).reshape(1, -1)
    B = np.exp(logB)
    positive = T > 0
    if not positive.any():
        raise ValueError('Empty injected template')
    scale = float(np.min(B[positive]/T[positive]))
    x = np.log10(scale)
    lo, hi = int(np.floor(x))-4, int(np.ceil(x))+3
    try:
        values = calc_max_lkl(Y, logB[None,:], pars, T[None,:],
                              np.ones(T.size), min_log10mu=lo,
                              max_log10mu=hi, n_seed_lkl=seeds)
        ll0, llmax, mu, sigma = map(jf, values)
        finite = all(v is not None for v in (ll0, llmax, mu, sigma))
        return {'finite': finite, 'mu_hat': mu, 'sigma_mu': sigma,
                'lnL0': ll0, 'lnLmax': llmax, 'seed_bounds': [lo, hi]}
    except (ValueError, ArithmeticError, AssertionError, FloatingPointError) as exc:
        return {'finite':False, 'error':repr(exc), 'seed_bounds':[lo,hi]}


def lkl_fixed(Y, logB, T, pars, mu):
    val = log_likelihood(mu, np.asarray(Y,float)[None,:], logB[None,:],
                         T[None,:], np.ones(len(T)), pars)
    if not np.isfinite(val):
        raise ValueError('Nonfinite original vendor likelihood')
    return float(val)


def analyze_one(case, channel, freq, indices, T, file, background, native_grid,
                trials, optimizer_pairs, seeds, rng, width, smoothing, chi_max,
                peak_over_background, halo_verification):
    with h5py.File(file, 'r') as h:
        required = {'frequency','PSD','logPSD','N_avg','N_kernel'}
        if not required.issubset(h.keys()):
            raise ValueError(f'{channel}: HDF5 missing native C columns')
        if len(h['frequency']) != native_grid.size:
            raise ValueError(f'{channel}: native C grid length mismatch')
        start = int(indices[0]); stop = int(indices[-1])+1
        if not np.array_equal(h['frequency'][start:stop], freq):
            raise ValueError(f'{channel}: exact native grid mismatch')
        measured_logP = np.asarray(h['logPSD'][start:stop],float)
        measured_P = np.asarray(h['PSD'][start:stop],float)
        n_avg = np.asarray(h['N_avg'][start:stop],float)
        if (not np.isfinite(measured_P).all() or np.any(measured_P<=0)
            or not np.isfinite(measured_logP).all() or np.any(n_avg<=0)):
            raise ValueError(f'{channel}: bad observed C spectrum')
        # Strictly check logPSD vs PSD, no reconstruction/interpolation.
        max_logdiff = float(np.max(np.abs(measured_logP-np.log(measured_P))))
        if max_logdiff > 1e-8:
            raise ValueError(f'{channel}: inconsistent measured PSD/logPSD: {max_logdiff}')
        center = float(freq[int(np.argmax(T))])
        lo = int(np.searchsorted(native_grid, center-width))
        hi = int(np.searchsorted(native_grid, center+width))
        regional = np.asarray(h['logPSD'][lo:hi],float)
    if regional.size <= max(100, 3*freq.size) or not np.isfinite(regional).all():
        raise ValueError('Too few valid neighboring bins to assess covariance')
    logB, B, pars, chi = vendor_fit(background, file, freq, chi_max)
    if not np.any(T>0) or np.any(T<0) or not np.isfinite(T).all():
        raise ValueError('Invalid signal template')
    ratio = float(np.max(T/B))
    if ratio <= 0 or not np.isfinite(ratio):
        raise ValueError('Signal-to-background ratio invalid')
    mu_inj = float(peak_over_background/ratio)
    # Retain correlations from actual data; remove slow trends, THEN use a
    # monotone quantile transformation to enforce fitted skew-normal marginal.
    hp = regional-gaussian_filter1d(regional, sigma=smoothing)
    rho = acf(hp, 20)
    ranks = rankdata(hp, method='average')
    uniform = np.clip((ranks - .5)/regional.size, 1e-10, 1-1e-10)
    loc, scale, alpha = [float(v) for v in pars[0]]
    mapped = skewnorm.ppf(uniform, a=alpha, loc=loc, scale=scale)
    if not np.isfinite(mapped).all():
        raise ValueError('Nonfinite skew-normal marginal transform')
    # Full blocks preserve the empirical 81-bin dependence; shuffled blocks
    # offer a reference null with the same marginal and NO adjacency structure.
    logLR_null, logLR_null_iid, logLR_inj = [], [], []
    optimizer = []
    for t in range(trials):
        start = int(rng.integers(0,len(mapped)-len(T)+1))
        block = mapped[start:start+len(T)]
        block_iid = mapped[rng.integers(0,len(mapped),size=len(T))]
        y_null = logB+block
        y_iid = logB+block_iid
        y_inj = np.log(np.exp(y_null)+mu_inj*T)
        # Call original vendor likelihood, not a replacement PDF/statistic.
        ln0 = lkl_fixed(y_null,logB,T,pars,0.0)
        ln1 = lkl_fixed(y_null,logB,T,pars,mu_inj)
        logLR_null.append(ln1-ln0)
        logLR_null_iid.append(lkl_fixed(y_iid,logB,T,pars,mu_inj)-
                              lkl_fixed(y_iid,logB,T,pars,0.0))
        logLR_inj.append(lkl_fixed(y_inj,logB,T,pars,mu_inj)-
                         lkl_fixed(y_inj,logB,T,pars,0.0))
        if t < optimizer_pairs:
            optimizer.append({'trial':t, 'null':fit_original(y_null, logB,T,pars,seeds),
                              'injected':fit_original(y_inj,logB,T,pars,seeds)})
    a,b,c = map(np.asarray, (logLR_null,logLR_null_iid,logLR_inj))
    fitted_inj = [o['injected'] for o in optimizer if o['injected']['finite']]
    go_signal = float(np.mean(c>0))
    go_null = float(np.mean(a<0))
    return {
        'case':case, 'channel':channel, 'frequency_hz':center,
        'n_bins':int(len(T)), 'actual_N_avg_peak':float(n_avg[np.argmax(T)]),
        'background_chi_square':chi, 'max_PSD_log_consistency_error':max_logdiff,
        'halo_qwp_recalculation_max_relative_error':halo_verification,
        'peak_T_over_B':ratio,
        'empirical_lag1':float(rho[1]), 'empirical_lag2':float(rho[2]),
        'linear_score_std_inflation_diagnostic':inflation_ratio(T,B,rho),
        'actual_neighborhood_bins':int(regional.size),
        'simulation':{'method':'contiguous-highpass-PSD empirical rank-copula, vendor skew-normal marginals',
                      'trials':trials, 'random_seed_deterministic':True,
                      'injected_mu':mu_inj,
                      'injected_peak_signal_to_background':peak_over_background,
                      'fixed_hypothesis_logLR':'original vendor logL(mu_injected)-logL(0)',
                      'null_median_logLR':float(np.median(a)),
                      'injected_median_logLR':float(np.median(c)),
                      'null_fraction_prefers_background':go_null,
                      'injected_fraction_prefers_signal':go_signal,
                      'null_logLR_std_correlated':float(np.std(a,ddof=1)),
                      'null_logLR_std_shuffled':float(np.std(b,ddof=1)),
                      'logLR_std_corr_to_shuffled':float(np.std(a,ddof=1)/np.std(b,ddof=1)) if np.std(b,ddof=1)>0 else None,
                      'note':'Fixed two-hypothesis LR is NOT profiled q0 or a p-value; marginal/covariance model is a surrogate.'},
        'original_vendor_optimizer_trials':optimizer,
        'optimizer_injections_finite':len(fitted_inj),
        'optimizer_injections_positive_fraction':float(np.mean([r['mu_hat']>0 for r in fitted_inj])) if fitted_inj else None,
        'signal_direction_check_pass':bool(go_signal>.75 and go_null>.75 and
                                           np.median(c)>np.median(a)),
        'vendor_optimizer_injection_check_pass':bool(optimizer_pairs == 0 or
                         any(row['mu_hat'] is not None and row['mu_hat'] > 0
                             for row in fitted_inj)),
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117'))
    p.add_argument('--qwp-npz',required=True,type=Path)
    p.add_argument('--halo-npz',required=True,type=Path)
    p.add_argument('--output-json',required=True,type=Path)
    p.add_argument('--trials',type=int,default=128)
    p.add_argument('--optimizer-pairs',type=int,default=2)
    p.add_argument('--seeds',type=int,default=21)
    p.add_argument('--seed',type=int,default=20261119)
    p.add_argument('--half-width-hz',type=float,default=20.)
    p.add_argument('--smooth-sigma-bins',type=float,default=100.)
    p.add_argument('--chi-max',type=float,default=10.)
    p.add_argument('--injected-peak-over-background',type=float,default=2.5)
    a=p.parse_args()
    if (a.trials < 32 or a.optimizer_pairs < 0 or a.optimizer_pairs>a.trials or
        a.seeds < 3 or a.half_width_hz<=0 or a.smooth_sigma_bins<=0 or
        a.chi_max<=0 or a.injected_peak_over_background <= 0):
        p.error('Invalid trial count or test controls')
    root=a.root.expanduser().resolve()
    out=a.output_json.expanduser().resolve()
    qwp=a.qwp_npz.expanduser().resolve()
    halo=a.halo_npz.expanduser().resolve()
    if out.suffix!='.json' or out.exists(): p.error('Output must be a NEW .json; no overwrites')
    if not qwp.is_file() or not halo.is_file(): p.error('QWP or halo diagnostic NPZ missing')
    spectra={c:root/'handoff'/f'APE_0_15435_{c}.h5' for c in ('X','Y')}
    backgrounds={c:root/'background'/f'background_{c}_FULL.json' for c in ('X','Y')}
    if any(not x.is_file() for x in [*spectra.values(),*backgrounds.values()]):
        p.error('Missing real X/Y HDF5 or vendor full fitted background')
    start_time=time.monotonic()
    rng=np.random.default_rng(a.seed)
    background={c:LPSDJSONIO(str(backgrounds[c])) for c in ('X','Y')}
    results=[]
    with np.load(qwp,allow_pickle=False) as z, np.load(halo,allow_pickle=False) as k:
        meta=json.loads(str(z['metadata_json'].item()))
        halo_meta=json.loads(str(k['metadata_json'].item()))
        if meta.get('status')!='PROVISIONAL_SELECTED_FREQUENCY_QWP_PSD_NOT_PHYSICAL_UPPER_LIMIT':
            p.error('Unexpected QWP template provenance')
        if halo_meta.get('status')!='DIAGNOSTIC_ONLY_NO_VALIDATED_AXION_LIMIT':
            p.error('Unexpected halo template provenance')
        if meta.get('source_halo_sha256')!=sha256(halo):
            p.error('QWP template is not derived from supplied halo file: SHA-256 mismatch')
        phys=meta.get('provisional_qwp_parameters',{})
        required={'finesse_s','tau_s','tau1_s','rho_dm_gev4','unit_conversion','reference_coupling_gev_inv'}
        if set(phys)!=required or any(not np.isfinite(float(v)) or float(v)<=0 for v in phys.values()):
            p.error('Incomplete provisional QWP working parameters')
        ids=sorted(int(key.split('_')[1]) for key in z.files
                   if key.startswith('case_') and key.endswith('_frequency_hz') and key.count('_')==3)
        if ids!=[0,1]:p.error('Expected exactly the existing 5 and 8 kHz cases')
        with h5py.File(spectra['X'],'r') as hx:
            native_grid=np.asarray(hx['frequency'][:],float)
            if (native_grid.size!=5035956 or not np.isfinite(native_grid).all()
                or np.any(np.diff(native_grid)<=0)):
                p.error('Unexpected production native-grid size or ordering')
        for case in ids:
            tag=f'case_{case}_'
            freq=np.asarray(z[tag+'frequency_hz'],float)
            idx=np.asarray(z[tag+'native_bin_indices'],np.int64)
            T=np.asarray(z[tag+'T_total_psd_per_unit_mu_rad2_hz'],float)
            halof=np.asarray(k[tag+'halo_frequency_hz'],float)
            hw=np.asarray(k[tag+'halo_weights'],float)
            kd=np.asarray(k[tag+'K_D_hz_inv_per_unit_variance'],float)
            if (len(freq)!=81 or idx.shape!=freq.shape or T.shape!=freq.shape
                or np.any(np.diff(idx)!=1) or np.any(np.diff(freq)<=0)
                or idx[0]<0 or idx[-1]>=native_grid.size
                or not np.array_equal(native_grid[idx[0]:idx[-1]+1],freq)
                or hw.shape!=halof.shape or kd.shape!=(81,len(halof))
                or not np.isclose(hw.sum(),1.,atol=1e-12,rtol=1e-12)):
                p.error(f'Case {case} failed native-grid/halo shape checks')
            for name in ('frequency_hz','native_bin_indices','halo_frequency_hz','halo_weights'):
                if not np.array_equal(z[tag+name],k[tag+name]):
                    p.error(f'Case {case}: QWP/halo input mismatch: {name}')
            R=qwp_amplitude_per_coupling(halof,finesse_s=float(phys['finesse_s']),
                     tau_s=float(phys['tau_s']),tau1_s=float(phys['tau1_s']),
                     rho_dm_gev4=float(phys['rho_dm_gev4']),
                     unit_conversion=float(phys['unit_conversion']))
            expected = .5*float(phys['reference_coupling_gev_inv'])**2 * (kd@(hw*R**2))
            rel=float(np.max(np.abs(T-expected))/np.max(expected))
            if rel>1e-12 or not np.isfinite(rel):
                p.error(f'Case {case}: saved QWP model fails independent template reconstruction')
            for ch in ('X','Y'):
                results.append(analyze_one(case,ch,freq,idx,T,spectra[ch],background[ch],
                                native_grid,a.trials,a.optimizer_pairs,a.seeds,rng,
                                a.half_width_hz,a.smooth_sigma_bins,a.chi_max,
                                a.injected_peak_over_background,rel))
                row=results[-1]
                print(f"case={case} {ch}: null prefers background "
                      f"{row['simulation']['null_fraction_prefers_background']:.1%}; "
                      f"injected prefers signal {row['simulation']['injected_fraction_prefers_signal']:.1%}; "
                      f"corr std~{row['linear_score_std_inflation_diagnostic']:.2f}",flush=True)
            center=float(freq[np.argmax(T)])
            low=int(np.searchsorted(native_grid,center-a.half_width_hz))
            high=int(np.searchsorted(native_grid,center+a.half_width_hz))
            local=[]
            for ch in ('X','Y'):
                with h5py.File(spectra[ch],'r') as hd:
                    q=np.asarray(hd['logPSD'][low:high],float)
                    local.append(q-gaussian_filter1d(q,sigma=a.smooth_sigma_bins))
            xy=float(np.corrcoef(local[0],local[1])[0,1])
            for row in results:
                if row['case']==case:row['empirical_XY_same_bin_correlation']=xy
            print(f"case={case}: observed X/Y residual correlation={xy:.4f}",flush=True)
    complete=all(r['signal_direction_check_pass'] and
                 r['vendor_optimizer_injection_check_pass'] for r in results)
    report={
        'status':'PASS_FAST_INTEGRATED_LIKELIHOOD_ENGINEERING' if complete else 'REVIEW_FAST_VALIDATION_DIAGNOSTICS',
        'science_status':'DIAGNOSTIC_ONLY_NOT_APPROVED_FOR_PHYSICAL_AXION_LIMITS',
        'go_for_more_diagnostic_likelihood_computations':bool(complete),
        'physical_inference_ready':False,
        'tested': ['saved halo->QWP->C D PSD formula','X/Y real native-grid and PSD/logPSD',
                   'original vendor fitted background','empirical 20Hz nearby-bin lag covariance',
                   'original vendor log-likelihood under null and power-domain injection',
                   'original signed-mu vendor optimizer on selected null/injection draws'],
        'provisional_assumptions':{
           'X_and_Y_each_alpha_equals_1':'SEPARATE COUNTERFACTUALS, NEVER JOINT',
           'halo':'benchmark fixed lab velocity and finite modes, not certified',
           'QWP_calibration':'working optical parameters; phase and X/Y projections unverified',
           'production_window_psll':'unverified',
           'Monte_Carlo':'empirical rank-copula based on one local real spectrum; not original C generated noise',
           'bin_correlations':'original vendor likelihood still assumes independence',
           'global_trials_and_95pct_coverage':'not established',
        },
        'inputs':{'qwp_npz':str(qwp),'qwp_sha256':sha256(qwp),
                  'halo_npz':str(halo),'halo_sha256':sha256(halo),
                  'spectra':{c:str(spectra[c]) for c in ('X','Y')},
                  'backgrounds':{c:str(backgrounds[c]) for c in ('X','Y')}},
        'controls':{'trials_per_case_channel':a.trials,'optimizer_pairs_per_case_channel':a.optimizer_pairs,
                    'vendor_seed_search':a.seeds,'random_seed':a.seed,
                    'injection_peak_to_background':a.injected_peak_over_background,
                    'half_width_hz':a.half_width_hz,
                    'smooth_sigma_bins':a.smooth_sigma_bins},
        'elapsed_seconds':time.monotonic()-start_time,
        'cases':results,
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x',encoding='utf-8') as fd:
        json.dump(report,fd,indent=2,allow_nan=False)
        fd.write('\n')
    print(json.dumps({'status':report['status'],'report':str(out),
       'elapsed_seconds':round(report['elapsed_seconds'],1),
       'summary':[{'case':r['case'],'channel':r['channel'],
                   'corr_std_ratio':round(r['linear_score_std_inflation_diagnostic'],3),
                   'null_correct':round(r['simulation']['null_fraction_prefers_background'],3),
                   'injected_correct':round(r['simulation']['injected_fraction_prefers_signal'],3),
                   'optimizer_injected_finite':r['optimizer_injections_finite'],
                   'X_Y_r':round(r['empirical_XY_same_bin_correlation'],4),
                   'signal_direction_check':r['signal_direction_check_pass'],
                   'vendor_optimizer_injection_check':r['vendor_optimizer_injection_check_pass']}
                  for r in results]},indent=2))


if __name__=='__main__':main()