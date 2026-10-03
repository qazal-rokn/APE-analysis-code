"""Path-only preflight tests; these intentionally never open a .dxd file."""
from pathlib import Path
import tempfile
import unittest

from ape_c_lpsd.ingestion.reference_paths import (
    ReferenceDewesoftPathError,
    resolve_reference_dewesoft_paths,
)


class ReferencePathResolutionTests(unittest.TestCase):
    def test_configured_reference_root_resolves_vendor_paths(self):
        root = Path("/Users/user/Desktop/APE Research/APE-data-analysis")
        resolved = resolve_reference_dewesoft_paths(root)
        self.assertEqual(resolved.reference_project_root, root.resolve())
        self.assertTrue(resolved.wrapper_path.is_file())
        self.assertTrue(resolved.library_path.is_file())

    def test_missing_path_error_identifies_root_and_wrapper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "missing-reference"
            with self.assertRaises(ReferenceDewesoftPathError) as caught:
                resolve_reference_dewesoft_paths(root)
        message = str(caught.exception)
        self.assertIn(str(root.resolve()), message)
        self.assertIn("DWDataReaderHeader.py", message)


if __name__ == "__main__":
    unittest.main()
