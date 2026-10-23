"""Deterministic calibrated-radian APE X/Y synthetic fixture.

This is a coherent-tone software smoke fixture, NOT a stochastic axion halo
injection. No calibration or physical-coupling recovery claims.
"""
import numpy as np
import h5py
from pathlib import Path

def make_fixture(destination, *, seed, frequency_hz, rms_rotation_rad,
                 noise_rms_rad, phase_rad=0., x_power_fraction=0.5,
                 fs_hz=2048., duration_s=60.):
    path=Path(destination).expanduser().resolve()
    if path.exists() or path.parent == path:
        raise FileExistsError("Refusing to overwrite synthetic input")
    if fs_hz!=2048. or duration_s!=60.:
        raise ValueError("Only 2048 Hz, 60-second smoke fixture is approved")
    vals=(frequency_hz,rms_rotation_rad,noise_rms_rad,phase_rad,x_power_fraction)
    if not all(np.isfinite(v) for v in vals) or not 5<=frequency_hz<=900 or rms_rotation_rad<0 or noise_rms_rad<=0 or not 0<=x_power_fraction<=1:
        raise ValueError("Invalid synthetic tone/noise parameters")
    rng=np.random.default_rng(int(seed))
    t=np.arange(int(fs_hz*duration_s))/fs_hz
    base=np.sqrt(2)*rms_rotation_rad*np.sin(2*np.pi*frequency_hz*t+phase_rad)
    x=base*np.sqrt(x_power_fraction)+rng.normal(0,noise_rms_rad,t.size)
    y=base*np.sqrt(1-x_power_fraction)+rng.normal(0,noise_rms_rad,t.size)
    path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,"x") as h5:
        h5.create_dataset("rho_x",data=x)
        h5.create_dataset("rho_y",data=y)
        h5.attrs["synthetic_only"]=True
        h5.attrs["unit"]="rad"
        h5.attrs["fs_hz"]=fs_hz
        h5.attrs["seed"]=int(seed)
        h5.attrs["coherent_tone_frequency_hz"]=float(frequency_hz)
        h5.attrs["coherent_tone_rms_rotation_rad"]=float(rms_rotation_rad)
    return dict(path=str(path),sample_count=len(t),signal_type="coherent_tone_not_axion_halo")
