"""Stochastic finite-duration axion halo polarization injection.

This generator accepts the *original* vendor local power shape weights and
frequencies. It does not invent an axion velocity distribution or a phase
transfer function. The QWP response is an amplitude envelope only; phase
response remains an explicit science-validation requirement.

Random independent complex Gaussian Fourier amplitudes generate a stationary
Gaussian realization with ensemble-mean power h_k times response_k**2.
"""
from pathlib import Path
import h5py
import numpy as np
from physics.qwp_signal import qwp_amplitude_per_coupling


def inject_discrete_halo(destination, *, fs_hz, duration_s, seed,
                         rest_frequency_hz, support_frequency_hz,
                         normalized_vendor_shape, coupling_gev_inv,
                         noise_rms_rad, quadrature_power_fraction,
                         finesse_s, tau_s, tau1_s, rho_dm_gev4,
                         unit_conversion):
    path=Path(destination).expanduser().resolve()
    if path.exists():
        raise FileExistsError("Refusing to overwrite synthetic halo fixture")
    fs=float(fs_hz); duration=float(duration_s)
    n=int(round(fs*duration))
    f=np.asarray(support_frequency_hz,dtype=float)
    weights=np.asarray(normalized_vendor_shape,dtype=float)
    a=np.asarray(quadrature_power_fraction,dtype=float)
    if (not np.isfinite(fs) or fs<=0 or not np.isfinite(duration) or duration<=0
            or n<4 or not np.isclose(n/fs,duration)
            or f.ndim!=1 or not len(f) or weights.shape!=f.shape
            or not np.isfinite(f).all() or np.any(f<=0)
            or np.any(f>=fs/2) or np.any(np.diff(f)<=0)
            or not np.isfinite(weights).all() or np.any(weights<0)
            or not np.isclose(weights.sum(),1.,rtol=1e-12,atol=1e-12)
            or a.shape!=(2,) or not np.isfinite(a).all()
            or np.any(a<0) or not np.isclose(a.sum(),1.)
            or not np.isfinite(coupling_gev_inv) or coupling_gev_inv<0
            or not np.isfinite(noise_rms_rad) or noise_rms_rad<=0):
        raise ValueError("Invalid discrete halo, sampling or X/Y inputs")
    grid=np.fft.rfftfreq(n,1/fs)
    idx=np.rint(f*duration).astype(int)
    if (np.any(idx<1) or np.any(idx>=len(grid)-1)
            or len(np.unique(idx))!=len(idx)
            or np.any(np.abs(grid[idx]-f)>1e-8*np.maximum(1.,f))):
        raise ValueError("Halo support must lie on distinct FFT bins; increase duration")
    if f[0]<rest_frequency_hz:
        raise ValueError("Signal support must be at/above axion rest frequency")
    rng=np.random.default_rng(int(seed))
    transfer=qwp_amplitude_per_coupling(
        f,finesse_s=finesse_s,tau_s=tau_s,tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4,unit_conversion=unit_conversion)
    # CN(0,1) amplitudes: E|z|^2=1.  For Fourier conventions,
    # a real sinusoid of RMS amplitude s has |rfft|=N*s/sqrt(2).
    z=(rng.normal(size=len(f))+1j*rng.normal(size=len(f)))/np.sqrt(2)
    base=np.zeros(len(grid),dtype=np.complex128)
    base[idx]=(n/np.sqrt(2))*float(coupling_gev_inv)*transfer*np.sqrt(weights)*z
    xcoeff=base*np.sqrt(a[0]); ycoeff=base*np.sqrt(a[1])
    x=np.fft.irfft(xcoeff,n=n)+rng.normal(0,noise_rms_rad,n)
    y=np.fft.irfft(ycoeff,n=n)+rng.normal(0,noise_rms_rad,n)
    path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,"x") as h5:
        h5.create_dataset("rho_x",data=x)
        h5.create_dataset("rho_y",data=y)
        h5.create_dataset("signal_frequency_hz",data=f)
        h5.create_dataset("vendor_shape_weights",data=weights)
        h5.attrs["synthetic_only"]=True
        h5.attrs["unit"]="rad"
        h5.attrs["model"]="stochastic_gaussian_vendor_shape"
        h5.attrs["seed"]=int(seed)
        h5.attrs["fs_hz"]=fs
        h5.attrs["duration_s"]=duration
        h5.attrs["rest_frequency_hz"]=float(rest_frequency_hz)
        h5.attrs["injected_coupling_gev_inv"]=float(coupling_gev_inv)
        h5.attrs["qwp_phase_response"]="not_included_envelope_only"
    return dict(path=str(path),n_samples=n,
                frequency_resolution_hz=1/duration,
                stochastic=True,phase_response_validated=False)
