#!/usr/bin/env python3
"""DIAGNOSTIC: local 81-bin provisional APE QWP+halo template in ORIGINAL vendor likelihood.

Uses previously generated provisional 5/8 kHz QWP NPZ and the real measured
4.28-hour X/Y PSD and fitted vendor background. Separately processes each
quadrature with a conditional UNIT projection alpha_channel = 1, not a joint
physical X/Y model.  No C reruns, no vendor-code modification, no physical UL.

IMPORTANT: Neighboring C Constant-Q bins overlap and are correlated.  The
unchanged vendor likelihood sums log-probabilities as if independent, so the
multi-bin likelihood results are a SOFTWARE / SHAPE-INTEGRATION DIAGNOSTIC,
NOT physically calibrated detection significance, uncertainty, or limits.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model
from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl


def finite_or_none(x):
    x = float(x)
    return x if np.isfinite(x) else None


def do_likelihood(root, channel, case, indices, f, T, vendor_background, seed_count, chi_cut):
    input_h5 = root / 'handoff' / f'APE_0_15435_{channel}.h5'
    if not input_h5.is_file():
        raise RuntimeError('Missing measured spectrum: '+str(input_h5))
    with h5py.File(input_h5, 'r') as file:
        n = len(file['frequency'])
        if indices[0] < 0 or indices[-1] >= n:
            raise RuntimeError(f'{channel} case {case}: native index out of range')
        first, last = int(indices[0]), int(indices[-1])+1
        if not np.array_equal(np.asarray(file['frequency'][first:last])[indices-first],f):
            raise RuntimeError(f'{channel} case {case}: native frequency grid mismatch')
        logP = np.asarray(file['logPSD'][first:last])[indices-first].astype(float)
        N = np.asarray(file['N_avg'][first:last])[indices-first].astype(float)
    if not np.isfinite(logP).all() or np.any(N<=0):
        raise RuntimeError(f'{channel} case {case}: invalid measured PSD')

    table = vendor_background.get_df(vendor_background.get_label(filepath=str(input_h5)))
    required={'fmin','fmax','x_knots','y_knots','alpha_skew','loc_skew','sigma_skew','chi_sqr'}
    if not required.issubset(table.columns):
        raise RuntimeError(f'{channel}: missing original vendor fitted-background columns')
    # Original vendor likelihood has one residual skew-normal parameter triple
    # per channel, NOT different parameters for bins straddling fit segments.
    eligible=table[(table['fmin'] <= float(f[0])) & (table['fmax'] >= float(f[-1]))]
    if eligible.empty:
        raise RuntimeError(f'{channel} case {case}: 81-bin region crosses original background-fit boundaries; cannot use one vendor residual model')
    fit=eligible.iloc[0]
    chi=float(fit['chi_sqr'])
    if not (np.isfinite(chi) and chi < chi_cut):
        raise RuntimeError(f'{channel} case {case}: background chi-square {chi} fails quality cut')
    knots=np.concatenate((np.asarray(fit['x_knots'],float),np.asarray(fit['y_knots'],float)))
    logB=np.asarray(bkg_model(np.log(f),knots),float).reshape(-1)
    if logB.shape != f.shape:
        raise RuntimeError(f'{channel} case {case}: incorrect spline output shape')
    B=np.exp(logB)
    pars=np.asarray([[fit['loc_skew'],fit['sigma_skew'],fit['alpha_skew']]],float)
    if not np.isfinite(B).all() or np.any(B<=0) or not np.isfinite(pars).all() or pars[0,1] <= 0:
        raise RuntimeError(f'{channel} case {case}: invalid fitted-background model')

    # All QWP/halo spectral power already resides in peak_norm = T_j.
    # Using peak_shape=ones prevents counting the halo twice. It ALSO makes
    # unchanged vendor get_valid_mu_ranges (which ignores peak_shape) valid.
    peak_norm=T[None,:]
    peak_shape=np.ones(f.size,float)
    # Scale the ORIGINAL optimizer's signed-mu search to the signal/background
    # ratio; never use LIGO's hard-coded mu ranges for APE.
    ratio=np.min(B/T[T>0]) if np.all(T>0) else np.min((B/T)[T>0])
    if not (np.isfinite(ratio) and ratio>0):
        raise RuntimeError(f'{channel} case {case}: nonfinite B/T scale')
    logscale=float(np.log10(ratio))
    lo=int(np.floor(logscale))-4
    hi=int(np.ceil(logscale))+3
    fit_result=calc_max_lkl(logP[None,:],logB[None,:],pars,
                            peak_norm,peak_shape,
                            min_log10mu=lo,max_log10mu=hi,n_seed_lkl=seed_count)
    ll0,llmax,muhat,sigma=[finite_or_none(z) for z in fit_result]
    ok=all(x is not None for x in (ll0,llmax,muhat,sigma))
    return {
      'case':case,'channel':channel,'native_bin_count':len(f),
      'frequency_interval_hz':[float(f[0]),float(f[-1])],
      'peak_native_frequency_hz':float(f[np.argmax(T)]),
      'max_template_to_background_PSD':float(np.max(T/B)),
      'fitted_background_chi_square':chi,
      'N_avg_at_peak':float(N[np.argmax(T)]),
      'projection_assumption':'alpha_this_channel=1 counterfactual; DO NOT combine with other channel',
      'vendor_peak_norm':'full native-bin QWP+halo T_j PSD per mu',
      'vendor_peak_shape':'ones; halo shape folded into peak_norm exactly once',
      'vendor_mu_seed_exponents':[lo,hi],
      'vendor_results':{'lnL0':ll0,'lnLmax':llmax,'mu_hat':muhat,'sigma_mu':sigma},
      'finite_numerical_result':ok,
      'no_physical_significance':True,
      'note':'Original vendor bins are treated independent, despite C overlap. This is solely a numerical interface exercise.',
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117'))
    p.add_argument('--qwp-npz',type=Path,required=True)
    p.add_argument('--output-json',type=Path,required=True)
    p.add_argument('--seeds',type=int,default=25)
    p.add_argument('--max-chi-sqr',type=float,default=10.)
    args=p.parse_args()
    root=args.root.expanduser().resolve()
    output=args.output_json.expanduser().resolve()
    if not args.qwp_npz.is_file():p.error('Missing QWP+halo NPZ')
    if output.suffix!='.json' or output.exists():p.error('Output must be new .json; no overwrites')
    if not (3<=args.seeds<=1000 and args.max_chi_sqr>0):p.error('Invalid seeds or chi-square cut')
    backgrounds={}
    for channel in ('X','Y'):
        bgfile=root/'background'/f'background_{channel}_FULL.json'
        if not bgfile.is_file():p.error(f'Missing background file: {bgfile}')
        backgrounds[channel]=LPSDJSONIO(str(bgfile))
    result=[]
    with np.load(args.qwp_npz,allow_pickle=False) as z:
        meta=json.loads(str(z['metadata_json'].item()))
        if meta.get('status')!='PROVISIONAL_SELECTED_FREQUENCY_QWP_PSD_NOT_PHYSICAL_UPPER_LIMIT':
            p.error('Input must be the known provisional QWP PSD template')
        ids=sorted({int(k.split('_')[1]) for k in z.files if k.startswith('case_') and k.endswith('_frequency_hz') and k.count('_')==3})
        if not ids:p.error('No saved QWP cases')
        for case in ids:
            prefix=f'case_{case}_'
            f=np.asarray(z[prefix+'frequency_hz'],float)
            indices=np.asarray(z[prefix+'native_bin_indices'],np.int64)
            T=np.asarray(z[prefix+'T_total_psd_per_unit_mu_rad2_hz'],float)
            if not (f.ndim==1 and f.size>=2 and indices.shape==T.shape==f.shape
                    and np.isfinite(f).all() and np.all(np.diff(f)>0)
                    and np.all(np.diff(indices)>0) and np.isfinite(T).all()
                    and np.all(T>=0) and T.max()>0):
                p.error(f'Bad QWP template case {case}')
            for channel in ('X','Y'):
                result.append(do_likelihood(root,channel,case,indices,f,T,
                                            backgrounds[channel],args.seeds,args.max_chi_sqr))
    status='PASS_NUMERICAL' if all(v['finite_numerical_result'] for v in result) else 'PARTIAL_NONFINITE_VENDOR_RESULT'
    out={
       'status':status+'_MULTIBIN_DIAGNOSTIC_ONLY',
       'input_qwp_npz':str(args.qwp_npz.expanduser().resolve()),
       'QWP_and_halo_physics':'PROVISIONAL_UNVERIFIED',
       'production_kaiser_psll':'UNVERIFIED',
       'X_Y_projection':'UNIT FOR EACH CHANNEL SEPARATELY, NOT SIMULTANEOUSLY',
       'statistical_validity':'INVALID FOR PHYSICAL INFERENCE: overlapping Constant-Q bins correlated; unchanged vendor likelihood assumes independence',
       'original_vendor_unchanged':True,
       'physical_g_upper_limit_computed':False,
       'diagnostic_cases':result,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as dst:
        json.dump(out,dst,indent=2,allow_nan=False)
        dst.write('\n')
    print(json.dumps({'status':out['status'],'cases':[
        {'case':a['case'],'channel':a['channel'],
         'bins':a['native_bin_count'],
         'max_T_over_B':a['max_template_to_background_PSD'],
         'mu_hat':a['vendor_results']['mu_hat'],
         'sigma_mu':a['vendor_results']['sigma_mu'],
         'finite':a['finite_numerical_result']}
        for a in result],'report':str(output)},indent=2))

if __name__=='__main__': main()