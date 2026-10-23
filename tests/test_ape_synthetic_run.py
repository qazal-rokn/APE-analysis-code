"""Tests of isolated synthetic fixture and PS batch join; no external C."""
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from ape_c_lpsd.analysis.synthetic_smoke import make_fixture
from tools.run_ape_synthetic_candidates import join_ps_sidecars


class Tests(unittest.TestCase):
    def test_synthetic_is_two_channels_and_immutable(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"fixture.h5"
            info=make_fixture(f,seed=1,frequency_hz=100.0,
                              rms_rotation_rad=1e-6,noise_rms_rad=1e-8)
            self.assertEqual(info["sample_count"],122880)
            with h5py.File(f) as h:
                self.assertTrue(h.attrs["synthetic_only"])
                self.assertEqual(h["rho_x"].shape,(122880,))
                self.assertEqual(h["rho_y"].shape,(122880,))
                self.assertFalse(np.array_equal(h["rho_x"][:],h["rho_y"][:]))
            with self.assertRaises(FileExistsError):
                make_fixture(f,seed=1,frequency_hz=100.0,
                             rms_rotation_rad=1e-6,noise_rms_rad=1e-8)

    def test_join_order_and_boundary_validation(self):
        with tempfile.TemporaryDirectory() as d:
            paths=[]
            for i,txt in enumerate(("1 2\n2 3\n","3 4\n4 5\n")):
                p=Path(d)/f"batch.{i}.ps"
                p.write_text(txt)
                paths.append(dict(start_index=i*2,rows=2,ps_file=str(p)))
            target=Path(d)/"combined.ps"
            self.assertEqual(join_ps_sidecars(paths,target),4)
            self.assertEqual(target.read_text(),"1 2\n2 3\n3 4\n4 5\n")
            with self.assertRaises(FileExistsError):
                join_ps_sidecars(paths,target)
            paths[1]["start_index"]=3
            with self.assertRaises(ValueError):
                join_ps_sidecars(paths,Path(d)/"wrong.ps")


if __name__=="__main__":
    unittest.main()
