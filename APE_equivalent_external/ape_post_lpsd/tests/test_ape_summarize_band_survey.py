"""Synthetic self-contained tests of survey consolidation and provenance."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ape_summarize_band_survey import collect,read_mass

class SurveySummaryTests(unittest.TestCase):
    def test_complete_and_partial_are_both_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)
            (p/"scan_manifest.json").write_text(json.dumps(dict(
                status="PROVISIONAL_SURVEY_NOT_PHYSICAL_LIMIT",
                config_hash="abc",rest_frequencies_hz=[100.,200.])))
            for i, status in enumerate(("complete","partial_sigma")):
                d=p/f"mass_{i:07d}";d.mkdir()
                obj=dict(index=i,rest_hz=[100.,200.][i],config_hash="abc",
                    X=dict(lnL0=-2.,lnLmax=-1.,mu_hat=2.,
                           sigma_mu=3.),
                    Y=dict(lnL0=-2.,lnLmax=-1.5,mu_hat=-1.,
                           sigma_mu=4. if i==0 else None))
                (d/("completed.json" if i==0 else "partial.json")).write_text(json.dumps(obj))
            _, rows,complete,partial=collect(p)
            self.assertEqual((len(rows),complete,partial),(4,1,1))
            self.assertEqual(rows[0]["q0_nominal_uncalibrated"],2.)
            self.assertEqual(rows[1]["q0_nominal_uncalibrated"],0.)
            self.assertIsNone(rows[-1]["sigma_mu_vendor_uncalibrated"])

    def test_provenance_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);(d/"completed.json").write_text(json.dumps(dict(index=0,
                rest_hz=100.,config_hash="wrong")))
            with self.assertRaises(ValueError):
                read_mass(d,0,100.,"expected")
