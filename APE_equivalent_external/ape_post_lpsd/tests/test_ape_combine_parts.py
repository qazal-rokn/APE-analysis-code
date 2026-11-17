"""Boundary tests for APE chunk provenance; vendor combiner isn't retested."""
import json
import tempfile
from pathlib import Path
from unittest import TestCase, mock
import unittest
import numpy as np

from ape_post_lpsd import ape_combine_parts as adapter
from ape_post_lpsd import ape_max_likelihood as likelihood


class TestAPEChunkIntegration(TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.root = Path(self.work.name)
        self.dir = self.root / 'parts'
        self.dir.mkdir()
        self.template = dict(channel='X', source_spectrum='APE_1_2_X.h5',
                             source_background='bkg.json', source_response='response.npz',
                             response_metadata={'mu_definition':'synthetic','response_units':'rad^2/Hz per unit mu'},
                             signal_model='single-bin additive PSD: exp(bkg)+mu*peak_norm; peak_shape=1',
                             frequency_range_hz=[65,8000], scan_log10mu=[-2,2],
                             joint_xy_likelihood_used=False, frequency_count=2)
        for n in range(2):
            path = self.dir / f'ape_X_{n}_result.npy'
            np.save(path, np.array([[65+n*2, -1., -0.8, .1, .2],
                                    [66+n*2, -1., -0.6, .2, .3]]))
            manifest = dict(self.template, start_index_within_analysis_range=n*2,
                            native_grid_indices=[n*2, (n+1)*2])
            path.with_suffix('.json').write_text(json.dumps(manifest))

    def test_delegates_once_to_vendor_combiner_and_checks_provenance(self):
        # This mock performs file assembly as a fixture; no testing of vendor
        # concatenation or sorting algorithm is needed here.
        def fake_vendor(*, results_dir, prefix, output_file, verbose):
            self.assertEqual(prefix, 'ape_X')
            self.assertFalse(verbose)
            arrays = [np.load(self.dir / f'ape_X_{n}_result.npy') for n in range(2)]
            np.save(output_file, np.concatenate(arrays))
        with mock.patch.object(adapter, 'vendor_combine', side_effect=fake_vendor) as vendor:
            manifest = adapter.combine_ape_parts(self.dir, 'ape_X', self.root/'all_X.npy')
        self.assertEqual(vendor.call_count, 1)
        self.assertEqual(manifest['frequency_count'], 4)
        self.assertEqual(manifest['native_grid_indices'], [0,4])
        self.assertEqual(len(manifest['combined_parts']), 2)
        self.assertTrue((self.root/'all_X.json').exists())

    def test_refuses_gapped_or_cross_channel_chunks(self):
        side = self.dir/'ape_X_1_result.json'
        d = json.loads(side.read_text())
        d['native_grid_indices'] = [5,7]
        side.write_text(json.dumps(d))
        with self.assertRaisesRegex(adapter.APECombineError, 'offset|gap'):
            adapter.combine_ape_parts(self.dir, 'ape_X', self.root/'bad.npy')
        self.assertFalse((self.root/'bad.npy').exists())
        d['native_grid_indices'] = [2,4]
        d['channel'] = 'Y'
        side.write_text(json.dumps(d))
        with self.assertRaisesRegex(adapter.APECombineError, 'Inconsistent APE'):
            adapter.combine_ape_parts(self.dir, 'ape_X', self.root/'bad.npy')

    def test_invalid_ape_chunk_offset_refused(self):
        side = self.dir/'ape_X_1_result.json'
        d = json.loads(side.read_text())
        d['start_index_within_analysis_range'] = 0
        side.write_text(json.dumps(d))
        with self.assertRaisesRegex(adapter.APECombineError, 'offset'):
            adapter.combine_ape_parts(self.dir, 'ape_X', self.root/'bad.npy')

    def test_stage4_new_only_chunk_option_available(self):
        self.assertIn('start_index', likelihood.run_single_quadrature.__code__.co_varnames)


if __name__=='__main__':
    unittest.main()
