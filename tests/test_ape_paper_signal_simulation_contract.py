"""Scientific guard: forbid stochastic DM signal in paper-based APE tests."""
import unittest
from ape_c_lpsd.analysis.paper_signal_simulation_contract import (
    validate_paper_signal_simulation_manifest as validate)

class TestPaperApproximation(unittest.TestCase):
    def setUp(self):
        self.valid = {
            "random_DM_phase": False,
            "random_DM_amplitude": False,
            "random_DM_velocity_modes": False,
            "stochasticity_of_DM_signal_neglected": True,
            "source_paper_doi": "10.1038/s41598-025-33428-2",
            "injection_domain": "AFTER_NATIVE_C_NOISE_PS_rad2_NOT_TIME_DOMAIN_DM",
        }
    def test_deterministic_power_accepted(self):
        self.assertTrue(validate(self.valid))
    def test_random_phase_refused(self):
        self.valid["random_DM_phase"] = True
        with self.assertRaises(ValueError):
            validate(self.valid)
    def test_random_amplitude_refused(self):
        self.valid["random_DM_amplitude"] = True
        with self.assertRaises(ValueError):
            validate(self.valid)
    def test_random_halo_mode_refused(self):
        self.valid["random_DM_velocity_modes"] = True
        with self.assertRaises(ValueError):
            validate(self.valid)
    def test_arbitrary_coherent_tone_is_not_validated_halo(self):
        self.valid["injection_domain"] = "coherent_250Hz_tone"
        with self.assertRaises(ValueError):
            validate(self.valid)
if __name__ == "__main__":
    unittest.main()
