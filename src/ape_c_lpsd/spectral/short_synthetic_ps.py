"""Isolated 60-second X/Y synthetic C-PS plan, separate from frozen APE profile.

Uses exactly the existing C _c_argv numerical argument path. Source data is
read-only, new output directories are exclusive, no real-data mode.
"""
from pathlib import Path
from ape_c_lpsd.spectral.constantq import _c_argv
from ape_c_lpsd.spectral.ps_optin import plan_ps_optin

def plan_short_synthetic(*, executable, ps_executable, ps_sha256,
                         input_hdf5, dataset, fft_hdf5, output_dir,
                         endpoints=(0,400,850,1301)):
    if dataset not in ("rho_x","rho_y"):
        raise ValueError("Synthetic dataset must be rho_x or rho_y")
    inp,fft,out=(Path(x).expanduser().resolve() for x in (input_hdf5,fft_hdf5,output_dir))
    if not inp.is_file() or not all(x.is_absolute() for x in (inp,fft,out)):
        raise ValueError("Existing absolute synthetic HDF5 input required")
    if not out.name.startswith("ape-synthetic-ps-") or out.exists():
        raise ValueError("Run output must be new, dedicated ape-synthetic-ps-* directory")
    if fft.parent != out:
        raise ValueError("FFT must be in new run output directory")
    if endpoints not in ((0,1301),(0,400,850,1301)):
        raise ValueError("Use exact separately approved 1301-bin test partition")
    args=dict(exe=executable,inp=inp,dataset="/"+dataset,fft=fft,
              fs=2048.,t_stop=60.,f_min=5.,f_max=900.,grid_bins=1301,
              overlap=20,regression=0,window=-2,psll=238.13,
              scale=1,epsilon=10,memory_power=26,
              kernel_threshold=1e-15,a_option=2,debug=False)
    fft_cmd=_c_argv(output=out/"fft_stage",method=3,start=0,count=1301,**args)
    cmds=[(i,_c_argv(output=out/f"{dataset}.lpsd.{i}",method=2,
                           start=start,count=stop-start,**args))
           for i,(start,stop) in enumerate(zip(endpoints,endpoints[1:]))]
    plan=dict(dataset=dataset,fft=fft_cmd,batches=cmds,output_dir=str(out),
              fft_hdf5=str(fft),c_config='GNUTERM "APE Synthetic Constant-Q"\nFORMAT "fDNbri"\nEND\n',
              format="fDNbri")
    return plan_ps_optin(plan,ps_executable=ps_executable,expected_sha256=ps_sha256)
