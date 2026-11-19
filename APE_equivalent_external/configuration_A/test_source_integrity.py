"""Regression guards for Configuration A. Run from APE_equivalent_external:
python -m unittest discover -s configuration_A -p 'test_*.py' -v
"""
import ast
import hashlib
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "scalar_dark_matter_lpsd"
EXPECTED = {
    "lvk/stats.py": "239ee4e4a00cff0a3bfe3aa16e22608d1b2d966c",
    "lvk/max_lkl_executable.py": "099a6c9b895bf8c0abb1a252a1f76c97667a08af",
}
FORBIDDEN = ("ape_axion_halo", "ape_native_kaiser", "ape_provisional_qwp_multibin", "configuration_B")


def git_blob_sha(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


class TestConfigurationAIsolation(unittest.TestCase):
    def test_original_vendor_scientific_files_unchanged(self):
        for relative_path, expected in EXPECTED.items():
            with self.subTest(relative_path=relative_path):
                self.assertEqual(git_blob_sha((VENDOR / relative_path).read_bytes()), expected)

    def test_configuration_a_does_not_import_multibin(self):
        for path in (ROOT / "configuration_A").glob("*.py"):
            with self.subTest(path=path):
                source = path.read_text(encoding="utf-8")
                ast.parse(source)
                for forbidden in FORBIDDEN:
                    self.assertNotIn(forbidden, source)

    def test_adapter_has_single_bin_shape(self):
        source = (ROOT / "ape_post_lpsd" / "ape_max_likelihood.py").read_text(encoding="utf-8")
        self.assertIn("peak_shape = np.ones(1, dtype=np.float64)", source)
        self.assertIn("from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl", source)


if __name__ == "__main__":
    unittest.main()
