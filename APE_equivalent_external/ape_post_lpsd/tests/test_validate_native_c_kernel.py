"""Synthetic tests for native-C D-column kernel validation (no C invoked)."""
import json
from pathlib import Path
import tempfile
import sys
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_native_c_kernel import (
    KernelValidationError, make_tone_fixture, validate_c_kernel, export_kernel, load_campaign,
)


class KernelTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.folder = Path(self.work.name)
        self.grid = np.array([4999.995, 5000., 5000.005, 5000.01])
        self.idx = np.array([0, 1, 3], dtype=np.int64)
        self.modes = np.array([5000., 5000.005])
        self.amps = [2.e-7, 5.e-7]
        self.expected = np.array([[[.1, .2], [8., 4.], [.6, .3]],
                                  [[.2, .4], [3., 2.], [1., .1]]])
        self.files = {'X': [], 'Y': []}
        for c, ch in enumerate(('X', 'Y')):
            for k in range(len(self.modes)):
                pair = []
                for amp_i, amp in enumerate(self.amps):
                    p = self.folder / f'{ch}_{k}_{amp_i}.lpsd'
                    result = np.zeros(len(self.grid))
                    result[self.idx] = self.expected[c, :, k] * (amp**2/2)
                    with p.open('x') as stream:
                        stream.write('# native C output format fDNbri\n')
                        for f, D in zip(self.grid, result):
                            stream.write(f'{f:.17g} {D:.17g} 2 11 0 0\n')
                    pair.append(str(p))
                self.files[ch].append(pair)
        self.options = dict(native_frequency_hz=self.grid, selected_indices=self.idx,
                            mode_frequencies_hz=self.modes, c_outputs=self.files,
                            peak_amplitudes_rad=self.amps, columns=6)

    def test_recovers_true_per_variance_psd_kernel_without_rescaling(self):
        kernel, diag = validate_c_kernel(**self.options)
        np.testing.assert_allclose(kernel, self.expected, rtol=2e-14)
        self.assertEqual(len(diag), 4)
        self.assertEqual(diag[0]['active_bins_tested'], 3)
        self.assertLess(diag[0]['max_fractional_kernel_disagreement'], 1e-14)

    def test_accepts_explicit_four_column_combined_c_outputs(self):
        copies = {'X': [], 'Y': []}
        for channel in ('X', 'Y'):
            for mode_i, paths in enumerate(self.files[channel]):
                pair = []
                for amp_i, old in enumerate(paths):
                    new = self.folder / f'combined_{channel}_{mode_i}_{amp_i}.lpsd'
                    rows = []
                    for line in Path(old).read_text().splitlines():
                        rows.append(line if line.startswith('#') else ' '.join(line.split()[:4]))
                    new.write_text('\n'.join(rows) + '\n')
                    pair.append(str(new))
                copies[channel].append(pair)
        kernel, _ = validate_c_kernel(**dict(self.options, c_outputs=copies, columns=4))
        np.testing.assert_allclose(kernel, self.expected, rtol=2e-14)

    def test_rejects_wrong_quadratic_scaling(self):
        p = Path(self.files['X'][0][1])
        lines = p.read_text().splitlines()
        parts = lines[2].split()
        parts[1] = str(float(parts[1]) * 1.1)
        lines[2] = ' '.join(parts)
        p.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(KernelValidationError, 'amplitude-squared parity'):
            validate_c_kernel(**self.options)

    def test_rejects_mismatched_native_frequencies(self):
        p = Path(self.files['Y'][0][0])
        p.write_text(p.read_text().replace('5000.005', '5000.006'))
        with self.assertRaisesRegex(KernelValidationError, 'frequency grid mismatch'):
            validate_c_kernel(**self.options)

    def test_rejects_metadata_mismatch(self):
        p = Path(self.files['Y'][0][0])
        p.write_text(p.read_text().replace(' 2 11 0 0', ' 3 11 0 0'))
        with self.assertRaisesRegex(KernelValidationError, 'averaging/kernel metadata'):
            validate_c_kernel(**self.options)

    def test_rejects_wrong_input_format(self):
        p = Path(self.files['X'][0][0])
        p.write_text(p.read_text().replace(' 0 0', ''))
        with self.assertRaisesRegex(KernelValidationError, 'expected 6 C columns'):
            validate_c_kernel(**self.options)

    def test_exports_step_five_ready_npz_with_provisional_status(self):
        out = self.folder / 'kernel.npz'
        manifest = export_kernel(out, **self.options,
                                 c_configuration_source='synthetic C flags fixture',
                                 c_executable_sha256='a'*64)
        with np.load(out, allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_frequency_hz'], self.grid)
            np.testing.assert_array_equal(f['output_bin_indices'], self.idx)
            np.testing.assert_allclose(f['native_c_psd_kernel_per_variance_hz'], self.expected)
            self.assertEqual(json.loads(str(f['metadata_json'].item())), manifest)
        self.assertIn('PROVISIONAL', manifest['status'])
        with self.assertRaisesRegex(KernelValidationError, 'overwrite'):
            export_kernel(out, **self.options, c_configuration_source='synthetic',
                          c_executable_sha256='a'*64)

    def test_campaign_reads_native_grid_from_existing_hdf5(self):
        source = self.folder / 'handoff.h5'
        with h5py.File(source, 'x') as h5:
            h5.create_dataset('frequency', data=self.grid)
        config = self.folder / 'campaign.json'
        config.write_text(json.dumps({'native_spectrum_hdf5': source.name,
                                      'selected_indices': self.idx.tolist()}))
        content = load_campaign(config)
        np.testing.assert_array_equal(content['native_frequency_hz'], self.grid)
        self.assertEqual(content['selected_indices'], self.idx.tolist())

    def test_real_tone_fixture_has_recorded_peak_amplitude_and_variance(self):
        file = self.folder / 'known_tone.h5'
        report = make_tone_fixture(file, frequency_hz=100., peak_amplitude_rad=4e-7,
                                   sample_rate_hz=2000., duration_s=0.1, channel='X')
        with h5py.File(file) as f:
            data = f['rho_x'][:]
            self.assertEqual(len(data), 200)
            np.testing.assert_allclose(data, 4e-7*np.cos(2*np.pi*100*np.arange(200)/2000), atol=1e-17)
            self.assertEqual(f.attrs['peak_amplitude_rad'], 4e-7)
            self.assertEqual(f.attrs['kind'], 'SYNTHETIC_TONE_ONLY_NOT_AXION')
        self.assertAlmostEqual(report['input_variance_rad2'], (4e-7)**2/2)
        with self.assertRaisesRegex(KernelValidationError, 'overwrite'):
            make_tone_fixture(file, frequency_hz=100., peak_amplitude_rad=4e-7,
                              sample_rate_hz=2000., duration_s=0.1, channel='X')

    def test_refuses_undocumented_binary_and_degenerate_grid(self):
        with self.assertRaisesRegex(KernelValidationError, 'SHA-256'):
            export_kernel(self.folder/'result.npz', **self.options,
                          c_configuration_source='synthetic', c_executable_sha256='not-a-hash')
        broken=dict(self.options, selected_indices=[0,0,3])
        with self.assertRaisesRegex(KernelValidationError, 'indices'):
            validate_c_kernel(**broken)

if __name__=='__main__':
    unittest.main()
