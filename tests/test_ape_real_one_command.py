"""Fast non-invasive guards for unified APE real-data orchestration."""
import json
import tempfile
import unittest
import hashlib
from pathlib import Path
from unittest.mock import patch

from tools import ape_real_one_command as integration


class TestRealOneCommand(unittest.TestCase):
    def test_matching_reused_native_ps(self):
        with tempfile.TemporaryDirectory() as tmp:
            spectra={}
            for name in ('rho_x','rho_y'):
                path=Path(tmp)/(name+'.ps')
                path.write_text('65.0 1e-12\n65.1 2e-12\n')
                spectra[name]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                               'unit':'rad^2','schema':'frequency_PS_2_columns'}
            cfg={'calibrated_hdf5_sha256':'f'*64,'spectra':spectra}
            self.assertEqual(set(integration.reuse_native(cfg,{'hdf5_sha256_recorded':'f'*64})),
                             {'rho_x','rho_y'})
            cfg['calibrated_hdf5_sha256']='0'*64
            with self.assertRaises(ValueError):
                integration.reuse_native(cfg,{'hdf5_sha256_recorded':'f'*64})

    def test_missing_c_settings_refused(self):
        with patch.object(integration,'docker_sha',return_value='f'*64):
            with self.assertRaisesRegex(ValueError,'settings'):
                integration.run_native_c(
                    {'container':'lpsd-linux','binary':'/work/lpsd-exec',
                     'binary_sha256':'f'*64,'c_settings':{}},
                    {'sample_rate_hz':20000,'samples_per_channel':308717000},
                    Path('/tmp/ape-nonexistent/native_c'),dry_run=True)

    def test_cannot_label_as_physical_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'config.json'
            path.write_text(json.dumps({'scientific_scope':'physical_published_exclusion'}))
            with self.assertRaisesRegex(ValueError,'not an authorized'):
                integration.run(path,dry_run=True)


if __name__ == '__main__':
    unittest.main()
