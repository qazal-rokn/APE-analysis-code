"""Minimal opt-in C PS adapter tests, no C execution or real data."""
import hashlib
import tempfile
import unittest
from pathlib import Path
from ape_c_lpsd.spectral.ps_optin import PSError, plan_ps_optin, split_ps_batch


class OptInPSTests(unittest.TestCase):
    def test_planner_preserves_original_and_hash_guard(self):
        with tempfile.TemporaryDirectory() as d:
            executable = Path(d)/"test-c"
            executable.write_bytes(b"synthetic fixture only")
            executable.chmod(0o700)
            command = [str(executable), "-o", str(Path(d)/"output.lpsd")]
            original = {"fft": command.copy(), "batches": [(0, command.copy())],
                        "c_config": 'FORMAT "fDNbri"\n', "format": "fDNbri"}
            digest = hashlib.sha256(executable.read_bytes()).hexdigest()
            result = plan_ps_optin(
                original, ps_executable=executable, expected_sha256=digest)
            self.assertEqual(original["format"], "fDNbri")
            self.assertEqual(original["batches"][0][1][-1], command[-1])
            self.assertEqual(result["format"], "fDSNbri")
            self.assertTrue(result["batches"][0][1][-1].endswith(".raw7"))
            with self.assertRaises(PSError):
                plan_ps_optin(original, ps_executable=executable,
                              expected_sha256="0"*64)

    def test_streaming_ps_legacy_columns_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            raw, old, ps = [Path(d)/name for name in ("raw", "legacy", "ps")]
            raw.write_text("1 2.00 3.00 4 5 6 7\n2 3 4 5 6 7 8\n")
            self.assertEqual(split_ps_batch(raw, old, ps, expected_rows=2), 2)
            self.assertEqual(old.read_text(),
                             "1 2.00 4 5 6 7\n2 3 5 6 7 8\n")
            self.assertEqual(ps.read_text(), "1 3.00\n2 4\n")
            with self.assertRaises(PSError):
                split_ps_batch(raw, old, ps, expected_rows=2)

    def test_invalid_frequency_rejects_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            raw, old, ps = [Path(d)/name for name in ("raw", "legacy", "ps")]
            raw.write_text("1 2 3 4 5 6 7\n1 2 3 4 5 6 7\n")
            with self.assertRaises(PSError):
                split_ps_batch(raw, old, ps, expected_rows=2)
            self.assertFalse(old.exists())
            self.assertFalse(ps.exists())


if __name__ == "__main__":
    unittest.main()
