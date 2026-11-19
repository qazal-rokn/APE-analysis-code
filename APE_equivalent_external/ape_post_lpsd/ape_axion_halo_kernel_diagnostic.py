#!/usr/bin/env python3
"""Diagnostic-only APE axion halo x analytic C-LPSD PSD kernel on native grid.

This does NOT represent a validated physical APE coupling response. Uses a
fixed-velocity lab-frame truncated Gaussian velocity benchmark, the previously
cross-checked analytic Kaiser estimator kernel and, only if explicitly supplied,
the existing provisional APE QWP Case-2 optical response.

No C jobs, no vendor edits, no background refits, no automatic correction to
make integrated C-LPSD PSD unity. Every target is an independent hypothetical
axion rest frequency (different axion mass), NOT two simultaneous modes.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import hashlib
from pathlib import Path
import numpy as np
import h5py

C = 299792458.0

def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fd:
        for data in iter(lambda: fd.read(1048576), b''):
            h.update(data)
    return h.hexdigest()

def import_kernel(p):
    spec = importlib.util.spec_from_file_location('ape_existing_native_kaiser_kernel', p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def halo_modes(*, rest, n_samples, n_modes, seed, sigma, escape, lab_speed):
    """Fixed-lab 3-D Gaussian velocity benchmark, truncated in GALACTIC frame.

    Maxwell SHM parameterization is not uniquely specified by sigma; must be
    reviewed for the physical analysis. This sampling is synthetic and is
    NOT a detector injection, and does not validate phase or covariance.
    """
    rng = np.random.default_rng(seed)
    observed = []
    required = n_samples
    while required:
        v = rng.normal(0., sigma, size=(max(2048, int(required*1.4)), 3))
        v = v[(v*v).sum(axis=1) <= escape*escape]
        if not len(v):
            raise ValueError('No halo velocities accepted')
        take = min(len(v), required)
        v = v[:take]
        v[:, 0] -= lab_speed
        observed.append(rest*(1.0+(v*v).sum(axis=1)/(2.0*C*C)))
        required -= take
    freq = np.concatenate(observed)
    # Bin on the actual physical sampled support, without discarding tails.
    edges = np.linspace(float(freq.min()), float(freq.max())+1e-10, n_modes+1)
    counts, edges = np.histogram(freq, edges)
    positive = counts > 0
    centers = (edges[:-1]+edges[1:])*0.5
    return centers[positive], counts[positive].astype(np.float64)/n_samples, {
        'sample_min_hz': float(freq.min()), 'sample_max_hz': float(freq.max()),
        'modes_nonzero': int(positive.sum()), 'mean_frequency_hz': float(freq.mean()),
        'weight_sum': float(np.sum(counts)/n_samples),
    }

def response_from_repo(frequencies, config):
    from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling
    required = {'finesse_s','tau_s','tau1_s','rho_dm_gev4','unit_conversion',
                'reference_coupling_gev_inv','quadrature_power_fraction',
                'calibration_version'}
    missing=required-set(config)
    if missing: raise ValueError(f'Missing explicit optical inputs: {sorted(missing)}')
    a = np.asarray(config['quadrature_power_fraction'],float)
    if a.shape!=(2,) or not np.isfinite(a).all() or (a<0).any() or (a>1).any() or not a.any():
        raise ValueError('Invalid X/Y quadrature fractions')
    if not isinstance(config['calibration_version'],str) or not config['calibration_version'].strip():
        raise ValueError('Missing calibration version')
    gref=float(config['reference_coupling_gev_inv'])
    if not np.isfinite(gref) or gref <= 0:raise ValueError('g_ref must be positive')
    R=qwp_amplitude_per_coupling(frequencies,**{k:config[k] for k in ('finesse_s','tau_s','tau1_s','rho_dm_gev4','unit_conversion')})
    return R,gref,a

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectrum', type=Path,required=True)
    p.add_argument('--native-kernel-py',type=Path,required=True,
                   help='Previously validated ape_native_kaiser_kernel.py; import it, do not replace it')
    p.add_argument('--rest-hz',type=float,action='append',required=True,
                   help='Independent hypothesized axion rest frequency; repeat as desired')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--physics-json',type=Path,
                   help='OPTIONAL: explicitly reviewed QWP params, gref, quadrature fractions. No defaults.')
    p.add_argument('--half-width-bins',type=int,default=40)
    p.add_argument('--halo-modes',type=int,default=64)
    p.add_argument('--halo-samples',type=int,default=120000)
    p.add_argument('--seed',type=int,default=20261119)
    p.add_argument('--sigma-m-s',type=float,default=229000*np.sqrt(3/2),
                   help='DIAGNOSTIC benchmark 1D galactic Gaussian velocity dispersion')
    p.add_argument('--escape-m-s',type=float,default=528000,
                   help='DIAGNOSTIC galactic escape speed; not laboratory escape speed')
    p.add_argument('--lab-speed-m-s',type=float,default=229000,
                   help='DIAGNOSTIC fixed laboratory speed; no annual modulation')
    p.add_argument('--psll',type=float,default=238.13,
                   help='PILOT Kaiser PSLL, still unverified for production')
    a=p.parse_args()
    if a.output.exists() or a.output.suffix!='.npz':p.error('Use a NEW .npz output, never overwrite')
    if not a.spectrum.is_file() or not a.native_kernel_py.is_file():p.error('Missing inputs')
    if not 1<=len(a.rest_hz)<=20:p.error('Supply 1 to 20 rest frequencies')
    if a.half_width_bins<5 or a.half_width_bins>5000:p.error('Invalid half-width')
    if a.halo_modes<8 or a.halo_modes>512 or a.halo_samples<1000:p.error('Invalid halo discretization')
    if any(x<=0 or not np.isfinite(x) for x in [a.sigma_m_s,a.escape_m_s,a.lab_speed_m_s]):p.error('Invalid velocities')
    if not a.sigma_m_s < C and a.escape_m_s < C and a.lab_speed_m_s < C:p.error('Superluminal velocities')
    kmod=import_kernel(a.native_kernel_py)
    physics=None
    if a.physics_json:
        physics=json.loads(a.physics_json.read_text())
        # Explicitly report this is still a provisional optical expression.
    datasets={}; summaries=[]
    with h5py.File(a.spectrum,'r') as h:
        ds=h['frequency']; attrs=h['logPSD'].attrs
        fs,fmin,fmax,jdes=(float(attrs['fs']),float(attrs['fmin']),float(attrs['fmax']),int(attrs['Jdes']))
        if len(ds)!=jdes:raise ValueError('Unexpected native grid size')
        for i,rest in enumerate(a.rest_hz):
            if not np.isfinite(rest) or rest<=0:raise ValueError('Invalid rest frequency')
            modes,w,halo_report=halo_modes(rest=rest,n_samples=a.halo_samples,n_modes=a.halo_modes,
                 seed=a.seed+i,sigma=a.sigma_m_s,escape=a.escape_m_s,lab_speed=a.lab_speed_m_s)
            j=int(round((jdes-1)*np.log(rest/fmin)/np.log(fmax/fmin)))
            if j<a.half_width_bins or j+a.half_width_bins>=jdes:raise ValueError('Rest frequency near edge of output grid')
            indices=np.arange(j-a.half_width_bins,j+a.half_width_bins+1,dtype=np.int64)
            f=np.asarray(ds[j-a.half_width_bins:j+a.half_width_bins+1],dtype=np.float64)
            if not np.allclose(f,fmin*np.exp(indices*np.log(fmax/fmin)/(jdes-1)),rtol=1e-8,atol=1e-7):
                raise ValueError('Native C frequency grid verification failed')
            K=kmod.model_kernel(f,modes,fs=fs,fmin=fmin,fmax=fmax,jdes=jdes,psll=a.psll)
            unit_variance=np.einsum('jk,k->j',K,w)
            if not np.isfinite(unit_variance).all() or not (unit_variance>=0).all():
                raise ValueError('Nonfinite/negative window convolution')
            prefix=f'case_{i}_'
            datasets[prefix+'frequency_hz']=f
            datasets[prefix+'native_bin_indices']=indices
            datasets[prefix+'halo_frequency_hz']=modes
            datasets[prefix+'halo_weights']=w
            datasets[prefix+'K_D_hz_inv_per_unit_variance']=K
            datasets[prefix+'halo_convolved_K_D_hz_inv']=unit_variance
            report={'rest_hz':rest,'center_native_hz':float(f[a.half_width_bins]),
                    'max_halo_convolved_psd_hz_inv':float(unit_variance.max()),
                    'peak_output_frequency_hz':float(f[np.argmax(unit_variance)]),**halo_report}
            if physics is not None:
                R,gref,alpha=response_from_repo(modes,physics)
                T=0.5*gref*gref * alpha[:,None] * np.einsum('jk,k->j',K,w*R*R)[None,:]
                if not np.isfinite(T).all() or (T<0).any() or not T.any():
                    raise ValueError('Optical response overflowed, underflowed or zero')
                datasets[prefix+'provisional_T_x_y_psd_per_mu_rad2_hz']=T
                report['provisional_T_X_max_rad2_hz_per_mu']=float(T[0].max())
                report['provisional_T_Y_max_rad2_hz_per_mu']=float(T[1].max())
            summaries.append(report)
    meta={'status':'DIAGNOSTIC_ONLY_NO_VALIDATED_AXION_LIMIT',
          'window':'analytic Kaiser PSD via previously crosschecked existing module',
          'spectrum_source':str(a.spectrum.resolve()),
          'kernel_source':str(a.native_kernel_py.resolve()),
          'kernel_sha256':sha256(a.native_kernel_py),
          'psll_db':a.psll,'production_psll_verified':False,
          'reference_optical_model':'ape_post_lpsd.physics.qwp_signal.qwp_amplitude_per_coupling',
          'optical_model_included':physics is not None,
          'physics_json_source':str(a.physics_json.resolve()) if a.physics_json else None,
          'velocity_benchmark':{'sigma_m_s':a.sigma_m_s,'escape_m_s':a.escape_m_s,
                                'lab_speed_m_s':a.lab_speed_m_s,'fixed_lab':True,
                                'galactic_cutoff':True,'annual_modulation':False},
          'halo_samples':a.halo_samples,'halo_modes':a.halo_modes,'seed':a.seed,
          'disclaimers':['Not a laboratory halo certification','No C frequency-bin covariance model',
                         'No C averaging covariance model','No measured X/Y optical phase or response',
                         'Not a production-full-band signal template','No artificial PSD integral renormalization']}
    datasets['metadata_json']=np.array(json.dumps(meta,sort_keys=True))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('xb') as fd:np.savez_compressed(fd,**datasets)
    print(json.dumps({'status':meta['status'],'physical_QWP_response_included':physics is not None,
          'cases':summaries,'output':str(a.output.resolve())},indent=2))

if __name__=='__main__':main()