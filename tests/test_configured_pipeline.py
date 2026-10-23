"""Synthetic-only controller safeguards, not a substitute for actual C execution."""
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from ape_c_lpsd.spectral.configured_ps import plan_configured_ps
from tools.run_ape_pipeline import validate_input, execute

class Tests(unittest.TestCase):
    def test_real_data_guard_refuses_before_any_write(self):
        with tempfile.TemporaryDirectory() as d:
            config=dict(synthetic_only=True,execute_c=True,
                        input_kind="calibrated_file",
                        input_hdf5=str(Path(d)/"nonexistent.h5"),
                        run_dir=str(Path(d)/"ape-pipeline-test"))
            with self.assertRaises((FileNotFoundError,ValueError)):
                execute(config)
            self.assertFalse(Path(config["run_dir"]).exists())

    def test_missing_synthetic_provenance_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"input.h5"
            with h5py.File(f,"x") as h:
                h.create_dataset("rho_x",data=np.ones(20))
                h.create_dataset("rho_y",data=np.ones(20))
                h.attrs["fs_hz"]=10.
                h.attrs["unit"]="rad"
            conf=dict(input_hdf5=str(f),input_kind="synthetic_halo_file",
                      c_settings={"fs_hz":10.,"duration_s":2.})
            with self.assertRaises(ValueError):
                validate_input(conf)

    def test_no_overwrite_of_run_directory(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"ape-pipeline-test"
            f.mkdir()
            with self.assertRaises(ValueError):
                execute(dict(synthetic_only=True,execute_c=True,
                             input_kind="calibrated_file",input_hdf5=str(f),
                             run_dir=str(f)))


if __name__=="__main__":
    unittest.main()
