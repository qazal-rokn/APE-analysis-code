"""Opt-in PS-enabled C plan for calibrated X/Y data of arbitrary duration.

This is a new run-specific plan. Frozen historical ConstantQProfile is untouched.
The operator supplies every numerical setting. No C approximation.
"""
import math
from pathlib import Path
from ape_c_lpsd.spectral.constantq import _c_argv
from ape_c_lpsd.spectral.ps_optin import plan_ps_optin

def plan_configured_ps(*, executable, ps_executable, ps_sha256,
                       input_hdf5, dataset, fft_hdf5, output_dir, settings):
    if dataset not in ("rho_x","rho_y"):
        raise ValueError("Only rho_x/rho_y calibrated radians supported")
    inp=Path(input_hdf5).expanduser().resolve()
    out=Path(output_dir).expanduser().resolve()
    fft=Path(fft_hdf5).expanduser().resolve()
    if not inp.is_file() or out.exists() or fft.parent!=out:
        raise ValueError("Input must exist; output must be new; FFT must be inside output")
    if not out.name.startswith("ape-ps-"):
        raise ValueError("Use a dedicated new ape-ps-* output directory")
    required=("fs_hz","duration_s","f_min_hz","f_max_hz","grid_bins",
              "overlap_percent","window","psll_db","regression",
              "input_scale","epsilon","memory_power","kernel_threshold",
              "a_option","batch_count")
    if set(settings)!=set(required):
        raise ValueError("All and only declared C numerical settings required: "+", ".join(required))
    s=settings
    fs,duration,fmin,fmax=(float(s[x]) for x in ("fs_hz","duration_s","f_min_hz","f_max_hz"))
    grid=int(s["grid_bins"]); batches=int(s["batch_count"])
    if (not all(math.isfinite(x) for x in (fs,duration,fmin,fmax))
        or not 0<fmin<fmax<fs/2 or duration<=0 or grid<2 or
        not 1<=batches<=grid or grid!=s["grid_bins"]
        or batches!=s["batch_count"]):
        raise ValueError("Invalid C frequency/grid/duration parameters")
    args=dict(exe=executable,inp=inp,dataset="/"+dataset,fft=fft,
              fs=fs,t_stop=duration,f_min=fmin,f_max=fmax,grid_bins=grid,
              overlap=int(s["overlap_percent"]),regression=int(s["regression"]),
              window=int(s["window"]),psll=float(s["psll_db"]),
              scale=float(s["input_scale"]),epsilon=int(s["epsilon"]),
              memory_power=int(s["memory_power"]),
              kernel_threshold=float(s["kernel_threshold"]),
              a_option=int(s["a_option"]),debug=False)
    # Pure frequency-index partition; never modifies C estimator's frequency grid.
    endpoints=[(i*grid)//batches for i in range(batches+1)]
    fftcmd=_c_argv(output=out/"fft_stage",method=3,start=0,count=grid,**args)
    batchcmd=[(i,_c_argv(output=out/f"{dataset}.lpsd.{i}",
                        method=2,start=endpoints[i],
                        count=endpoints[i+1]-endpoints[i],**args))
              for i in range(batches)]
    legacy=dict(dataset=dataset,fft=fftcmd,batches=batchcmd,
                output_dir=str(out),fft_hdf5=str(fft),
                c_config='GNUTERM "APE configured PS"\nFORMAT "fDNbri"\nEND\n',
                format="fDNbri")
    return plan_ps_optin(legacy,ps_executable=ps_executable,
                         expected_sha256=ps_sha256)
