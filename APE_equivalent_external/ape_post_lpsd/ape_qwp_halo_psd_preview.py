#!/usr/bin/env python3
"""Provisional APE QWP + halo + native C-LPSD PSD kernel at selected masses.

Read the already-generated 5/8 kHz diagnostic halo kernel, apply EXACTLY the
existing APE qwp_amplitude_per_coupling function and compare the expected
signal PSD to real X/Y PSD samples at the same native C-LPSD bins.

NO full-band coverage, no calibration certification, no X/Y signal fractions,
no likelihood or physical coupling limits. Input PSD/halo files are read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--halo-npz', required=True, type=Path)
    p.add_argument('--spectrum-x', required=True, type=Path)
    p.add_argument('--spectrum-y', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    # These previously established APE *working* values are provisional.
    # They must not be silently equated with measured optical calibration.
    p.add_argument('--finesse-s', type=float, default=1500.0)
    p.add_argument('--tau-s', type=float, default=2.5e-9)
    p.add_argument('--tau1-s', type=float, default=1e-8)
    p.add_argument('--rho-dm-gev4', type=float, default=3.1e-42)
    p.add_argument('--unit-conversion', type=float, default=1.519e24)
    p.add_argument('--g-ref-gev-inv', type=float, default=1e-10)
    args = p.parse_args()

    inputs = (args.halo_npz, args.spectrum_x, args.spectrum_y)
    for source in inputs:
        if not source.is_file():
            p.error(f'Missing file: {source}')
    if args.output.suffix != '.npz' or args.output.exists():
        p.error('Output must be a new .npz; refusing to overwrite')
    phys = dict(finesse_s=args.finesse_s, tau_s=args.tau_s,
                tau1_s=args.tau1_s, rho_dm_gev4=args.rho_dm_gev4,
                unit_conversion=args.unit_conversion,
                reference_coupling_gev_inv=args.g_ref_gev_inv)
    if any(not np.isfinite(v) or v <= 0 for v in phys.values()):
        p.error('Every working physics parameter must be finite and positive')

    outputs = {}
    reports = []
    with (np.load(args.halo_npz, allow_pickle=False) as z,
          h5py.File(args.spectrum_x, 'r') as x,
          h5py.File(args.spectrum_y, 'r') as y):
        source_metadata = json.loads(str(z['metadata_json'].item()))
        if source_metadata.get('optical_model_included'):
            p.error('Kernel input already contains optical response; refusing to apply it twice')
        if source_metadata.get('status') != 'DIAGNOSTIC_ONLY_NO_VALIDATED_AXION_LIMIT':
            p.error('Unexpected source kernel status')
        cases = []
        for key in z.files:
            if (key.startswith('case_') and key.endswith('_frequency_hz')
                    and key.count('_') == 3):
                cases.append(int(key.split('_')[1]))
        cases.sort()
        if not cases or cases != list(range(len(cases))):
            p.error('Malformed or missing case_0..N kernel arrays')

        for channel, h in (('X',x),('Y',y)):
            if not all(field in h for field in ('PSD','frequency','N_avg')):
                p.error(f'{channel}: missing required HDF5 datasets')
            if h['PSD'].shape != h['frequency'].shape:
                p.error(f'{channel}: inconsistent PSD/frequency sizes')
            if len(h['PSD']) != len(x['PSD']):
                p.error('X/Y native grids have different sizes')
        for number in cases:
            tag=f'case_{number}_'
            f=np.asarray(z[tag+'frequency_hz'], dtype=np.float64)
            idx=np.asarray(z[tag+'native_bin_indices'], dtype=np.int64)
            halo=np.asarray(z[tag+'halo_frequency_hz'], dtype=np.float64)
            weights=np.asarray(z[tag+'halo_weights'], dtype=np.float64)
            k=np.asarray(z[tag+'K_D_hz_inv_per_unit_variance'], dtype=np.float64)
            folded=np.asarray(z[tag+'halo_convolved_K_D_hz_inv'], dtype=np.float64)
            if not (idx.ndim == 1 and f.shape == idx.shape and len(idx)>0
                    and np.all(np.diff(idx)>0) and idx[0]>=0
                    and idx[-1] < len(x['PSD'])):
                p.error(f'{tag}: invalid native bin indices')
            if not (halo.ndim==1 and weights.shape==halo.shape
                    and k.shape==(len(f),len(halo)) and folded.shape==f.shape):
                p.error(f'{tag}: incompatible input kernel shapes')
            if (not np.isfinite(k).all() or (k < 0).any()
                    or not np.isfinite(weights).all() or (weights < 0).any()
                    or not np.isclose(weights.sum(),1.,rtol=1e-12,atol=1e-12)):
                p.error(f'{tag}: invalid PSD kernel or halo weights')
            recalculated=k@weights
            if not np.allclose(recalculated,folded,rtol=2e-10,atol=1e-12):
                p.error(f'{tag}: saved halo-convolved kernel is inconsistent')
            for channel,h in (('X',x),('Y',y)):
                native=h['frequency'][int(idx[0]):int(idx[-1])+1]
                if not np.array_equal(native, f):
                    p.error(f'{tag}: frequency grid mismatch in channel {channel}')
                PSD=np.asarray(h['PSD'][int(idx[0]):int(idx[-1])+1],float)
                N_avg=np.asarray(h['N_avg'][int(idx[0]):int(idx[-1])+1],float)
                if not np.isfinite(PSD).all() or (PSD <= 0).any():
                    p.error(f'{channel}: nonpositive PSD around {f[len(f)//2]} Hz')
                outputs[tag+f'observed_PSD_{channel}_rad2_hz']=PSD
                outputs[tag+f'N_avg_{channel}']=N_avg
            R=qwp_amplitude_per_coupling(halo, finesse_s=args.finesse_s,
                    tau_s=args.tau_s, tau1_s=args.tau1_s,
                    rho_dm_gev4=args.rho_dm_gev4,
                    unit_conversion=args.unit_conversion)
            # unit-variance PSD response K_D * [mean squared amplitude / 2]
            T=0.5*(args.g_ref_gev_inv**2) * (k @ (weights * R**2))
            if not np.isfinite(T).all() or (T < 0).any() or T.max() <= 0:
                p.error(f'{tag}: invalid optical PSD template')
            peak=int(np.argmax(T))
            outputs[tag+'frequency_hz']=f
            outputs[tag+'native_bin_indices']=idx
            outputs[tag+'halo_frequency_hz']=halo
            outputs[tag+'halo_weights']=weights
            outputs[tag+'R_QWP_rad_per_GeV_inv']=R
            outputs[tag+'T_total_psd_per_unit_mu_rad2_hz']=T
            # IMPORTANT: T_total corresponds to a unit quadrature fraction.
            # It is NOT valid to assign this full power to both X and Y.
            reports.append({
              'case':number,
              'center_native_frequency_hz':float(f[len(f)//2]),
              'peak_output_hz':float(f[peak]),
              'peak_total_PSD_rad2_hz_at_gref_for_alpha_equal_1':float(T[peak]),
              'observed_PSD_X_rad2_hz_at_peak':float(outputs[tag+'observed_PSD_X_rad2_hz'][peak]),
              'observed_PSD_Y_rad2_hz_at_peak':float(outputs[tag+'observed_PSD_Y_rad2_hz'][peak]),
              'signal_to_observed_PSD_X_if_alpha_equal_1':float(T[peak]/outputs[tag+'observed_PSD_X_rad2_hz'][peak]),
              'signal_to_observed_PSD_Y_if_alpha_equal_1':float(T[peak]/outputs[tag+'observed_PSD_Y_rad2_hz'][peak]),
              'N_avg_X_at_peak':float(outputs[tag+'N_avg_X'][peak]),
              'N_avg_Y_at_peak':float(outputs[tag+'N_avg_Y'][peak]),
            })

    meta = {
      'status': 'PROVISIONAL_SELECTED_FREQUENCY_QWP_PSD_NOT_PHYSICAL_UPPER_LIMIT',
      'provisional_qwp_parameters':phys,
      'mu_definition':'(g_a_gamma_gamma/g_ref)^2',
      'total_template_definition':'T_total = (g_ref^2 / 2) * sum_k halo_weight_k * R_QWP(f_k)^2 * K_D(output_bin,mode_k)',
      'quadrature_definition':'T_X = alpha_X*T_total and T_Y = alpha_Y*T_total; alpha_X/alpha_Y unknown and NOT assigned',
      'source_halo_npz':str(args.halo_npz.resolve()),
      'source_halo_sha256':digest(args.halo_npz),
      'x_spectrum':str(args.spectrum_x.resolve()),
      'y_spectrum':str(args.spectrum_y.resolve()),
      'QWP_code':'original ape_post_lpsd.physics.qwp_signal.qwp_amplitude_per_coupling',
      'production_Kaiser_PSLL_certified':False,
      'halo_is_provisional':True,
      'overlapping_bin_correlation_model_certified':False,
      'notes':'Only a 5/8 kHz local-band diagnostic. Observed PSD is not background fit or observed upper limit. No X/Y signal projection is inferred.'
    }
    outputs['metadata_json']=np.array(json.dumps(meta,sort_keys=True))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('xb') as fd:
        np.savez_compressed(fd,**outputs)
    print(json.dumps({'status':meta['status'],
           'physics_working_values':phys, 'unknown_X_Y_fractions':True,
           'cases':reports,'output':str(args.output.resolve())},indent=2))

if __name__=='__main__':
    main()