import json, tempfile
from pathlib import Path
import numpy as np
from ape_c_lpsd.core import export_hdf5, longest_segment, derived_grid, preflight_nodes, ValidationError

def test_hdf5_order_and_units():
    with tempfile.TemporaryDirectory() as d:
        x=np.arange(17,dtype=float); assert export_hdf5(x,Path(d)/"a.h5","/strain",20.,"rho_x","rad")

def test_longest_cross_block_without_gap_stitching():
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)
        blocks=[]
        for i,m in enumerate(([1,1,0,1],[1,1,1,0])):
            np.save(p/f"b{i}.npy",np.array(m,dtype=bool)); np.save(p/f"x{i}.npy",np.arange(4)); np.save(p/f"y{i}.npy",np.arange(4))
            blocks.append({"block_id":i,"start_idx":4*i,"stop_idx":4*i+4,"valid_file":f"b{i}.npy","rho_x_file":f"x{i}.npy","rho_y_file":f"y{i}.npy"})
        (p/"cache_metadata.json").write_text(json.dumps({"sample_rate_hz":2,"blocks":blocks}))
        _,s=longest_segment(p); assert (s.start,s.stop,s.count)==(3,7,4)

def test_grid_is_metadata_derived():
    fmin,fmax,j=derived_grid(1000,20000,1e6); assert fmin==20 and fmax<10000 and j>0

def test_preflight_refuses_infeasible_node():
    try: preflight_nodes(128,128,9,36,16,{"window":-2})
    except ValidationError as e: assert "required_dft_window_samples" in str(e)
    else: raise AssertionError("infeasible C grid was accepted")
