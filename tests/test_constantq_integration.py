"""Small, no-vendor/no-real-data checks of historical command and batch contracts."""
import tempfile
import unittest
from pathlib import Path
from ape_c_lpsd.spectral.constantq import (
    ConstantQError, ConstantQProfile, batches, plan_commands, FREQUENCIES
)
from ape_c_lpsd.spectral.combine import validate_and_combine, check_xy_frequency_alignment
from ape_c_lpsd.backends.slurm import render_slurm_preview


class ControllerTests(unittest.TestCase):
    def test_exact_historical_partition(self):
        groups = batches()
        self.assertEqual(len(groups), 16)
        self.assertEqual(sum(n for _, _, n in groups), FREQUENCIES)
        self.assertEqual(groups[0], (0, 0, 3649978))

    def test_command_parity_and_fft_once(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            p = plan_commands(root/"lpsd-exec", root/"input.h5", "rho_x",
                              root/"fft.h5", root/"new-run")
            self.assertEqual(len(p["batches"]), 16)
            self.assertEqual(p["fft"][p["fft"].index("-h")+1], "3")
            for _, argv in p["batches"]:
                self.assertEqual(argv[argv.index("-h")+1], "2")
                self.assertEqual(argv[argv.index("-O")+1], str(root/"fft.h5"))
            render = render_slurm_preview(p)
            self.assertIn("--dependency=afterok", render["required_dependency"])
            self.assertIn("-h 3", render["fft"])
            self.assertIn("-h 2", render["batches"][0])
            self.assertLess(render["fft"].index("#SBATCH"), render["fft"].index("set -euo"))
            self.assertLess(render["batches"][0].index("#SBATCH"), render["batches"][0].index("set -euo"))

    def test_refuses_magnitude_and_changed_profile(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            with self.assertRaises(ConstantQError):
                plan_commands(root/"c", root/"in.h5", "rho_mag", root/"fft.h5", root/"run")
            with self.assertRaises(ConstantQError):
                ConstantQProfile(overlap=30).validate_reference()

    def test_combine_and_reject_missing_bins(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            a, b = root/"0", root/"1"
            a.write_text("1 4 2 3 0 0\n2 5 2 3 0 0\n")
            b.write_text("3 6 3 4 0 0\n")
            out = root/"combined"
            self.assertEqual(validate_and_combine([a,b],out,endpoints=(0,2,3)), 3)
            self.assertEqual(len(out.read_text().splitlines()), 3)
            with self.assertRaises(ConstantQError):
                validate_and_combine([a,b],root/"bad",endpoints=(0,2,4))
            self.assertFalse((root/"bad").exists())

    def test_xy_rejects_trailing_row(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            x, y = root/"x", root/"y"
            x.write_text("1 2 3 4\n2 3 4 5\n")
            y.write_text("1 2 3 4\n2 3 4 5\n")
            self.assertEqual(check_xy_frequency_alignment(x,y), 2)
            y.write_text("1 2 3 4\n2 3 4 5\n3 4 5 6\n")
            with self.assertRaises(ConstantQError):
                check_xy_frequency_alignment(x,y)


if __name__ == "__main__":
    unittest.main()
