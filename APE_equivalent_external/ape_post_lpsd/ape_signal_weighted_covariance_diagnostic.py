#!/usr/bin/env python3
"""Non-invasive covariance diagnostic for previously computed APE QWP+halo templates.

Reads the observed native-grid PSDs, original fitted vendor backgrounds and the
provisional 81-bin axion signal T_j. Estimates empirical lag correlations of
high-pass logPSD residuals across +/- 20 Hz, and evaluates a *linearized*,
score-weighted variance inflation relative to independent native bins.

This is NOT the covariance of the original vendor skew-normal score, NOT an
uncertainty correction to be applied to vendor output, NOT coverage, and NOT
an axion upper limit. No original vendor estimator, optimizer, or fit changes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d

from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model


def lag_corr(z: np.ndarray, maxlag: int) -> np.ndarray:
    """Detrended logPSD Pearson ACF; no row/column covariance matrix allocated."""
    x = np.asarray(z, dtype=float).reshape(-1)
    if x.size < maxlag * 20 or not np.isfinite(x).all():
        raise ValueError('Insufficient finite logPSD samples for autocorrelation')
    x = x - np.mean(x)
    d = float(np.dot(x, x))
    if not np.isfinite(d) or d <= 0:
        raise ValueError('Zero/nonfinite logPSD residual variance')
    # A standard lag-correlation diagnostic on the same detrended residuals as
    # the prior command. The denominator is recomputed for each shifted pair.
    v = np.ones(maxlag+1)
    for k in range(1, maxlag+1):
        a,b = x[:-k],x[k:]
        aa=a-a.mean();bb=b-b.mean()
        den=np.linalg.norm(aa)*np.linalg.norm(bb)
        v[k]=float(np.dot(aa,bb)/den) if den else 0.
    return v


def weighted_diagnostic(w: np.ndarray, rho: np.ndarray) -> dict:
    """Variance ratio for an illustrative linear statistic with weights T/B.

    Numerator = w^T R w, R_ij=rho[|i-j|], zero for unmeasured lags.
    Denominator = w^T w (independent, equal residual marginal variance).
    This is not the exact vendor skew-normal likelihood's information matrix.
    """
    w=np.asarray(w,float)
    if w.ndim != 1 or w.size < 3 or not np.isfinite(w).all() or np.any(w < 0) or not np.any(w):
        raise ValueError('Need finite nonnegative signal/background weights')
    iid=float(np.dot(w,w))
    cov=iid
    for k in range(1,min(len(rho),len(w))):
        cov += 2.*float(rho[k])*float(np.dot(w[:-k],w[k:]))
    if cov <= 0 or not np.isfinite(cov):
        raise ValueError('Estimated lag covariance gives nonpositive weighted variance')
    n_eff_iid=float(w.sum()**2/iid)
    n_eff_corr=float(w.sum()**2/cov)
    return {'linear_score_variance_ratio_correlated_to_iid':cov/iid,
            'linear_score_std_ratio_correlated_to_iid':float(np.sqrt(cov/iid)),
            'signal_weighted_effective_bins_if_iid':n_eff_iid,
            'signal_weighted_effective_bins_with_empirical_lags':n_eff_corr,
            'score_weighting':'w_j = T_j / B_j; illustrative linearized logPSD score, NOT vendor skewnormal score'}


def fit_bg(bg: LPSDJSONIO, spectrum: Path, f: np.ndarray) -> np.ndarray:
    tab=bg.get_df(bg.get_label(filepath=str(spectrum)))
    required={'fmin','fmax','x_knots','y_knots','chi_sqr'}
    if not required.issubset(tab.columns):
        raise ValueError('Original background record lacks required columns')
    valid=tab[(tab['fmin']<=f[0]) & (tab['fmax']>=f[-1])]
    if valid.empty:
        raise ValueError('Case crosses fitted-background segment boundary')
    fit=valid.iloc[0]
    if not (np.isfinite(fit['chi_sqr']) and float(fit['chi_sqr'])<10.):
        raise ValueError('Background fit does not meet chi-square cut')
    knots=np.concatenate((np.asarray(fit['x_knots'],float),np.asarray(fit['y_knots'],float)))
    b=np.exp(np.asarray(bkg_model(np.log(f),knots),float).reshape(-1))
    if b.shape!=f.shape or not np.isfinite(b).all() or np.any(b<=0):
        raise ValueError('Invalid fitted background spectrum')
    return b


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117'))
    p.add_argument('--qwp-npz',required=True,type=Path)
    p.add_argument('--output-json',required=True,type=Path)
    p.add_argument('--half-width-hz',type=float,default=20.)
    p.add_argument('--smooth-sigma-bins',type=float,default=100.)
    p.add_argument('--max-lag',type=int,default=20)
    a=p.parse_args()
    if a.half_width_hz<=0 or a.smooth_sigma_bins<=0 or not 1<=a.max_lag<=50:
        p.error('Invalid diagnostic settings')
    root=a.root.expanduser().resolve()
    out=a.output_json.expanduser().resolve()
    if out.suffix!='.json' or out.exists():p.error('Output must be new .json (no overwrite)')
    if not a.qwp_npz.is_file():p.error('QWP template missing')
    h5={c:root/'handoff'/f'APE_0_15435_{c}.h5' for c in ('X','Y')}
    bg={c:root/'background'/f'background_{c}_FULL.json' for c in ('X','Y')}
    if any(not v.is_file() for v in (*h5.values(),*bg.values())):
        p.error('Missing original X/Y handoff HDF5 or background_FULL.json')
    backgrounds={c:LPSDJSONIO(str(bg[c])) for c in ('X','Y')}
    result=[]
    with np.load(a.qwp_npz,allow_pickle=False) as z:
        meta=json.loads(str(z['metadata_json'].item()))
        if meta.get('status')!='PROVISIONAL_SELECTED_FREQUENCY_QWP_PSD_NOT_PHYSICAL_UPPER_LIMIT':
            p.error('Unexpected QWP template provenance')
        ids=sorted(int(key.split('_')[1]) for key in z.files
                   if key.startswith('case_') and key.endswith('_frequency_hz') and key.count('_')==3)
        if not ids:p.error('No 81-bin cases saved in QWP template')
        for case in ids:
            pre=f'case_{case}_'
            f=np.asarray(z[pre+'frequency_hz'],float)
            indices=np.asarray(z[pre+'native_bin_indices'],np.int64)
            T=np.asarray(z[pre+'T_total_psd_per_unit_mu_rad2_hz'],float)
            if (f.ndim!=1 or len(f)<3 or indices.shape!=f.shape or T.shape!=f.shape
                or np.any(np.diff(indices)!=1) or not np.all(np.diff(f)>0)
                or not np.isfinite(T).all() or not np.any(T>0) or np.any(T<0)):
                raise ValueError(f'Case {case} bad native grid/template')
            center=float(f[np.argmax(T)])
            residuals={}
            rhos={}
            for ch in ('X','Y'):
                with h5py.File(h5[ch],'r') as h:
                    native=h['frequency']
                    lo=int(np.searchsorted(native[:],center-a.half_width_hz))
                    hi=int(np.searchsorted(native[:],center+a.half_width_hz))
                    if not np.array_equal(np.asarray(native[indices[0]:indices[-1]+1]),f):
                        raise ValueError(f'{ch} case {case}: native frequency mismatch')
                    PSD=np.asarray(h['PSD'][lo:hi],float)
                    if np.any(PSD<=0) or not np.isfinite(PSD).all():
                        raise ValueError(f'{ch} case {case}: invalid PSD in neighborhood')
                    logPSD=np.log(PSD)
                if len(logPSD)<a.max_lag*20:
                    raise ValueError('Not enough native bins for covariance diagnostic')
                resid=logPSD-gaussian_filter1d(logPSD,sigma=a.smooth_sigma_bins)
                residuals[ch]=resid
                rho=lag_corr(resid,a.max_lag)
                rhos[ch]=rho
                B=fit_bg(backgrounds[ch],h5[ch],f)
                d=weighted_diagnostic(T/B,rho)
                result.append({
                    'case':case,'channel':ch,'template_peak_frequency_hz':center,
                    'neighborhood_bins':len(logPSD),'lag_correlations':{
                        str(k):float(rho[k]) for k in [1,2,3,4,5,10,20] if k<len(rho)},
                    'max_T_over_B':float(np.max(T/B)),
                    **d,
                    'interpretation':'PROVISIONAL; not a corrected sigma_mu, significance or physical UL'})
            xy=float(np.corrcoef(residuals['X'],residuals['Y'])[0,1])
            for entry in result:
                if entry['case']==case:
                    entry['xy_residual_correlation_same_bin']=xy
    report={
        'status':'EMPIRICAL_CORRELATED_BIN_DIAGNOSTIC_NOT_VENDOR_LIKELIHOOD_CALIBRATION',
        'data_source':'Existing 4.28h APE X/Y PSDs plus original fitted vendor backgrounds',
        'signal':'Saved provisional 5/8kHz QWP+halo T_j, unit projection into each channel separately',
        'lag_method':'Gaussian-detrended log(PSD) in +/-20 Hz neighborhood; stationary Toeplitz approximation',
        'lag_truncation_warning':'Lags beyond --max-lag treated as zero for the illustrative weighted variance; no exact covariance/coverage implied',
        'vendor_likelihood_modified':False,
        'physics_inference_authorized':False,
        'cases':result}
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x',encoding='utf-8') as fh:json.dump(report,fh,indent=2,allow_nan=False);fh.write('\n')
    print(json.dumps({'status':report['status'],
       'cases':[{'case':r['case'],'channel':r['channel'],
                 'lag1':r['lag_correlations']['1'],'lag2':r['lag_correlations']['2'],
                 'xy_r':r['xy_residual_correlation_same_bin'],
                 'linearized_std_ratio':r['linear_score_std_ratio_correlated_to_iid'],
                 'effective_bins_iid':r['signal_weighted_effective_bins_if_iid'],
                 'effective_bins_correlated':r['signal_weighted_effective_bins_with_empirical_lags']}
                for r in result], 'report':str(out)},indent=2))

if __name__=='__main__':main()