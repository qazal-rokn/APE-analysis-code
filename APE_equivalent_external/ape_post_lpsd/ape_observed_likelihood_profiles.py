#!/usr/bin/env python3
"""Observed APE QWP/halo likelihood profiles at SAVED 5/8 kHz cases.

Directly evaluates the UNMODIFIED vendor stats.log_likelihood at each trial
mu = (g/g_ref)^2, for separate X and Y measured spectra.  Reads original
vendor fitted background and previously computed original vendor optimizer
results; does NOT refit background or rerun C-LPSD or recompute axion template.

IMPORTANT: This is a conditional/composite diagnostic, NOT a calibrated
physical upper limit or discovery significance.  C-LPSD bins are correlated,
X/Y signal fractions unknown, halo/QWP/production PSLL provisional.  Negative
mu values are included only to view the *signed* vendor fit shape.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.stats import log_likelihood
from ape_onepass_likelihood_validation import vendor_fit, fit_original


def original_logl(mu: float, observed_logpsd, log_background, signal, residual_params):
    """Call vendor implementation, never reimplement its PDF or optimization."""
    return float(log_likelihood(
        float(mu), observed_logpsd[None, :], log_background[None, :],
        signal[None, :], np.ones(signal.size), residual_params
    ))


def load_existing_fits(path: Path):
    if not path.is_file():
        return {}
    with path.open() as source:
        obj=json.load(source)
    if not obj.get('status', '').startswith('PASS_NUMERICAL'):
        raise ValueError('Existing optimizer results are not marked PASS_NUMERICAL')
    return {(int(x['case']), x['channel']): x['vendor_results']
            for x in obj['diagnostic_cases']}


def calculate(root: Path, template: Path, outdir: Path, existing_results: Path,
              seeds: int, chi_max: float, points: int, make_plot: bool):
    fits = load_existing_fits(existing_results)
    bgs={q:LPSDJSONIO(str(root/'background'/f'background_{q}_FULL.json'))
         for q in ('X','Y')}
    results=[]
    plot_rows=[]
    with np.load(template, allow_pickle=False) as z:
        metadata=json.loads(str(z['metadata_json'].item()))
        if metadata.get('status') != 'PROVISIONAL_SELECTED_FREQUENCY_QWP_PSD_NOT_PHYSICAL_UPPER_LIMIT':
            raise ValueError('Input is not the saved provisional APE QWP PSD')
        gref=float(metadata['provisional_qwp_parameters']['reference_coupling_gev_inv'])
        case_ids=sorted(int(k.split('_')[1]) for k in z.files
                 if k.startswith('case_') and k.endswith('_frequency_hz') and k.count('_')==3)
        if not case_ids:
            raise ValueError('No signal cases available')
        for case in case_ids:
            prefix=f'case_{case}_'
            f=np.asarray(z[prefix+'frequency_hz'],dtype=float)
            idx=np.asarray(z[prefix+'native_bin_indices'],dtype=np.int64)
            T=np.asarray(z[prefix+'T_total_psd_per_unit_mu_rad2_hz'],dtype=float)
            if (f.ndim!=1 or len(f)<2 or np.any(np.diff(f)<=0)
                or idx.shape!=f.shape or np.any(np.diff(idx)!=1)
                or T.shape!=f.shape or not np.isfinite(T).all()
                or np.any(T<0) or not np.any(T>0)):
                raise ValueError(f'Case {case} template failed validation')
            for q in ('X','Y'):
                spectrum=root/'handoff'/f'APE_0_15435_{q}.h5'
                with h5py.File(spectrum, 'r') as h:
                    start,stop=int(idx[0]),int(idx[-1])+1
                    actual_f=np.asarray(h['frequency'][start:stop],float)
                    Y=np.asarray(h['logPSD'][start:stop],float)
                    N=np.asarray(h['N_avg'][start:stop],float)
                if (not np.array_equal(actual_f,f) or not np.isfinite(Y).all()
                    or np.any(N<=0)):
                    raise ValueError(f'{case}/{q}: measured PSD/native grid invalid')
                logB,B,pars,chi=vendor_fit(bgs[q], spectrum, f, chi_max)
                rmax=float(np.max(T/B))
                scale=1/rmax    # s=mu/scale is peak signal/background ratio.
                if not np.isfinite(scale) or scale<=0:
                    raise ValueError('Invalid mu scale')
                key=(case,q)
                fit=fits.get(key)
                if fit is None:
                    fit=fit_original(Y,logB,T,pars,seeds)
                    if not fit['finite']:
                        raise ValueError(f'{case}/{q}: original vendor optimization failed: {fit}')
                muh=float(fit['mu_hat'])
                sig=float(fit['sigma_mu'])
                if not np.isfinite(muh) or not np.isfinite(sig) or sig <= 0:
                    raise ValueError(f'{case}/{q}: invalid saved vendor fit')
                # Evaluate a dense interval around the fitted best point and the
                # physical mu=0 boundary, keeping B+mu*T>0 for all native bins.
                best_s=muh/scale
                std_s=sig/scale
                right=max(2.,best_s+5*abs(std_s),5*abs(std_s))
                right=min(right,100.)
                left=max(-.98,min(-.2,best_s-5*abs(std_s)))
                pos=np.unique(np.r_[np.linspace(0.,min(2.,right),points//2),
                                      np.linspace(min(2.,right),right,points//4),
                                      best_s if best_s>=0 else 0.])
                neg=np.linspace(left,0.,points//4,endpoint=False)
                s=np.sort(np.unique(np.r_[neg,pos,0.]))
                mus=s*scale
                logl=np.array([original_logl(mu,Y,logB,T,pars) for mu in mus])
                if not np.isfinite(logl).all():
                    raise ValueError(f'{case}/{q}: vendor likelihood nonfinite')
                lkl0=original_logl(0.,Y,logB,T,pars)
                llmax=float(fit['lnLmax'])
                if abs(lkl0-float(fit['lnL0']))>1e-5:
                    raise ValueError(f'{case}/{q}: saved vendor lnL0 inconsistent with live likelihood')
                optimizer_profile_disagreement=float(max(0.,logl.max()-llmax))
                if optimizer_profile_disagreement>0.01:
                    print(f'WARNING {case}/{q}: cached vendor optimizer maximum is below sampled profile by '
                          f'{optimizer_profile_disagreement:.4g} log-likelihood units',flush=True)
                # The nonnegative-mu profile is diagnostic only. q0 does NOT
                # follow a known chi-square null law with these correlated bins.
                positive_max=max(lkl0,float(np.max(logl[mus>=0])))
                nominal_q0=max(0.,2*(positive_max-lkl0))
                maxref=max(llmax,float(logl.max()))
                arr=np.column_stack([mus,s,logl,2*(maxref-logl)])
                csv_path=outdir/f'case_{case}_{q}_vendor_likelihood_profile.csv'
                np.savetxt(csv_path,arr,delimiter=',',comments='',
                           header='mu,peak_signal_over_background,original_vendor_lnL,minus2DeltaLnL_signed_fit',
                           fmt='%.16g')
                row={
                  'case':case, 'channel':q,
                  'rest_frequency_approx_hz':float(np.round(f[len(f)//2])),
                  'signal_peak_frequency_hz':float(f[np.argmax(T)]),
                  'native_bins':len(f), 'g_ref_gev_inv':gref,
                  'fitted_background_chi_square':float(chi),
                  'N_avg_at_signal_peak':float(N[np.argmax(T)]),
                  'peak_T_over_B_per_mu':rmax,
                  'mu_hat_original_vendor_signed':muh,
                  'sigma_mu_original_vendor_uncalibrated':sig,
                  'lnL0_original_vendor':lkl0,
                  'lnLmax_original_vendor_signed':llmax,
                  'cached_optimizer_profile_disagreement_lnL':optimizer_profile_disagreement,
                  'nominal_q0_nonnegative_mu_UNCALIBRATED':nominal_q0,
                  'profile_csv':str(csv_path),
                  'likelihood_interpretation':'diagnostic original vendor composite likelihood; correlated bins, no calibrated p-value/UL',
                  'projection_assumption':'alpha_channel=1 separately for X and Y; NOT joint'
                }
                results.append(row)
                plot_rows.append((case,q,mus,logl,maxref,muh))
                print(f"case {case} {q}: mu_hat={muh:.6g}, sigma_vendor={sig:.6g}, "
                      f"nominal_q0_uncalibrated={nominal_q0:.4g}; profile saved",flush=True)
    plot_path=None
    if make_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig,axes=plt.subplots(len(case_ids),2,figsize=(11,3.6*len(case_ids)),squeeze=False)
            for (case,q,mus,logl,maxref,muh) in plot_rows:
                ix=case_ids.index(case)
                ax=axes[ix,0 if q=='X' else 1]
                ax.plot(mus/1e6,2*(maxref-logl),lw=1.5)
                ax.axvline(0.,ls='--',lw=1,color='gray')
                ax.axvline(muh/1e6,ls=':',lw=1.1,color='gray')
                ax.set_ylim(0,20)
                ax.set_xlabel(r'$\mu$ ($10^6$ units)')
                ax.set_ylabel(r'$-2\Delta\ln\mathcal{L}$ (signed-fit ref.)')
                peak_hz=next(r['signal_peak_frequency_hz'] for r in results if r['case']==case and r['channel']==q)
                ax.set_title(f'{peak_hz:.3f} Hz | {q} | observed vendor likelihood')
                ax.grid(alpha=.2)
            fig.suptitle('APE 5/8 kHz: original vendor likelihood (provisional, correlated bins)',fontsize=12)
            fig.tight_layout()
            plot_path=outdir/'ape_observed_vendor_likelihood_profiles.png'
            fig.savefig(plot_path,dpi=160)
            plt.close(fig)
        except (ImportError,OSError,RuntimeError) as exc:
            print('Plot skipped (JSON and CSV complete):',repr(exc),flush=True)
    report={
      'status':'OBSERVED_ORIGINAL_VENDOR_LIKELIHOOD_PROFILES_COMPUTED',
      'science_status':'PROVISIONAL_DIAGNOSTIC_NOT_CALIBRATED_FOR_PHYSICAL_SIGNIFICANCE_OR_LIMITS',
      'note':'This uses the original vendor likelihood at real observed PSD points; confidence/coverage NOT established.',
      'original_vendor_changed':False,
      'input_qwp':str(template),
      'cached_original_vendor_fit_file':str(existing_results) if fits else None,
      'optical_and_halo_model':'provisional, X/Y fraction unmeasured',
      'frequency_bin_covariance':'NONZERO; UNCHANGED VENDOR LIKELIHOOD ASSUMES INDEPENDENCE',
      'profiles':results,
      'plot':str(plot_path) if plot_path is not None else None,
    }
    with (outdir/'observed_likelihood_summary.json').open('x') as fh:
        json.dump(report,fh,indent=2,allow_nan=False)
        fh.write('\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117'))
    p.add_argument('--qwp-npz',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--existing-results-json',type=Path,default=None)
    p.add_argument('--seeds',type=int,default=25)
    p.add_argument('--max-chi-square',type=float,default=10.)
    p.add_argument('--points',type=int,default=400)
    p.add_argument('--no-plot',action='store_true')
    a=p.parse_args()
    root=a.root.expanduser().resolve()
    source=a.qwp_npz.expanduser().resolve()
    outdir=a.output_dir.expanduser().resolve()
    existing=(a.existing_results_json or root/'diagnostics'/'ape_qwp_vendor_multibin_5k_8k.json').expanduser().resolve()
    if not source.is_file() or a.points<100 or a.seeds<3:
        p.error('Missing QWP file or invalid numerical controls')
    if outdir.exists():p.error('Output directory already exists; refusing overwrite')
    outdir.mkdir(parents=True)
    result=calculate(root,source,outdir,existing,a.seeds,a.max_chi_square,a.points,not a.no_plot)
    print(json.dumps({'status':result['status'],
         'summary_json':str(outdir/'observed_likelihood_summary.json'),
         'plot_png':result['plot'],
         'profiles':[{k:r[k] for k in ('case','channel','mu_hat_original_vendor_signed','sigma_mu_original_vendor_uncalibrated','nominal_q0_nonnegative_mu_UNCALIBRATED')} for r in result['profiles']]},indent=2))


if __name__=='__main__':main()