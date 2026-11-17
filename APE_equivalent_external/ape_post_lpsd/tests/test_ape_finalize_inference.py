"""Tests ONLY the new APE boundary; vendor scientific functions stay unchanged."""
import json
import tempfile
from pathlib import Path
from unittest import TestCase, mock
import unittest

import numpy as np

from ape_post_lpsd import ape_finalize_inference as ape


class APEDownstreamInterface(TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.data = np.array([
            [100., -5., -5., -0.01, .2],
            [101., -5., 5., 0.2, .3],
            [102., -5., 11., 0.4, .4],
        ])
        self.npy = self.root / 'X_maxlkl.npy'
        np.save(self.npy, self.data)
        self.meta = {'response_units': 'rad^2/Hz per unit mu',
                     'mu_definition': 'mu=(g_a_gamma_gamma/reference_coupling_gev_inv)^2',
                     'reference_coupling_gev_inv': 1e-10,
                     'response_source': 'synthetic', 'optical_response_status': 'PROVISIONAL'}
        self.sidecar = self.root / 'X_maxlkl.json'
        self.save_manifest()

    def save_manifest(self):
        self.sidecar.write_text(json.dumps({
            'channel': 'X', 'joint_xy_likelihood_used': False,
            'signal_model': 'single-bin additive PSD: exp(bkg)+mu*peak_norm; peak_shape=1',
            'response_metadata': self.meta}))

    def test_vendor_screening_retains_trailing_cluster(self):
        result, info = ape.original_candidate_screen(self.data, sigma_limit=1.)
        self.assertEqual(result.shape, (1, 4))
        self.assertEqual(result[0, 0], 102.)
        self.assertEqual(result[0, 2], .4)
        self.assertEqual(info['clusters_retained'], 1)
        self.assertFalse(info['cross_quadrature_veto_applied'])

    def test_uncertainty_forwarded_to_original_asymptotic_function(self):
        with mock.patch.object(ape, 'get_upper_limits_approx', return_value=(np.array([.7]), np.array([.08]))) as original:
            output = ape.original_expected_ul(self.data, g_ref=1e-10, alpha_cl=.9, tol=1e-6)
        self.assertEqual(original.call_count, 3)
        original.assert_any_call(.2, alpha_CL=.9, tol=1e-6)
        np.testing.assert_allclose(output[:, 1], .7)
        np.testing.assert_allclose(output[:, 3], 7e-11)
        np.testing.assert_allclose(output[:, 4], 8e-12)

    def test_skips_invalid_sigma_without_changing_vendor_solver(self):
        data = self.data.copy()
        data[1, 4] = np.nan
        with mock.patch.object(ape, 'get_upper_limits_approx', return_value=(np.array([.5]), np.array([.1]))) as original:
            output = ape.original_expected_ul(data)
        self.assertEqual(original.call_count, 2)
        self.assertTrue(np.isnan(output[1, 1]))
        self.assertTrue(np.isnan(output[:, 3]).all())

    def test_complete_ape_software_outputs_are_provisional(self):
        out = self.root / 'X_diagnostic'
        with mock.patch.object(ape, 'get_upper_limits_approx', return_value=(np.array([.5]), np.array([.1]))):
            manifest = ape.finalize_single_quadrature(self.npy, self.sidecar, out, sigma_limit=1.)
        self.assertFalse(manifest['physics_validated'])
        self.assertFalse(manifest['observed_upper_limit_calculated'])
        cand = np.loadtxt(str(out) + '_screened_candidates.csv', delimiter=',', comments='#', ndmin=2)
        lim = np.loadtxt(str(out) + '_expected_ul.csv', delimiter=',', comments='#', ndmin=2)
        self.assertEqual(cand.shape, (1, 4))
        self.assertEqual(lim.shape, (3, 5))
        self.assertEqual(lim[0, 3], 5e-11)
        self.assertEqual(json.loads(Path(str(out) + '_analysis_status.json').read_text()), manifest)
        with self.assertRaisesRegex(ape.APEInferenceError, 'overwrite'):
            ape.finalize_single_quadrature(self.npy, self.sidecar, out)

    def test_no_coupling_units_without_matching_metadata(self):
        self.meta['mu_definition'] = 'synthetic arbitrary mu'
        self.save_manifest()
        _, _, gref = ape._read_validated_inputs(self.npy, self.sidecar)
        self.assertIsNone(gref)

    def test_joint_quadrature_and_unsorted_rows_rejected(self):
        record = json.loads(self.sidecar.read_text())
        record['joint_xy_likelihood_used'] = True
        self.sidecar.write_text(json.dumps(record))
        with self.assertRaisesRegex(ape.APEInferenceError, 'Joint'):
            ape._read_validated_inputs(self.npy, self.sidecar)
        record['joint_xy_likelihood_used'] = False
        self.sidecar.write_text(json.dumps(record))
        arr = self.data[::-1]
        np.save(self.npy, arr)
        with self.assertRaisesRegex(ape.APEInferenceError, 'ordered'):
            ape._read_validated_inputs(self.npy, self.sidecar)

    def test_no_implicit_LIGO_consistency_or_strain_conversion(self):
        self.assertNotIn('get_A_star', ape.__dict__)
        self.assertNotIn('analyze_per_ifo', ape.__dict__)
        self.assertNotIn('get_beta', ape.__dict__)
        self.assertNotIn('RHO_LOCAL', ape.__dict__)


if __name__ == '__main__':
    unittest.main()
