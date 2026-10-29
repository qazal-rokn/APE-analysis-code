"""Fail-closed 2026-paper deterministic axion-power simulation contract.

This is validation only; it changes no C, halo, QWP, or vendor algorithm.
"""
from __future__ import annotations

REQUIRED_FIXED_FLAGS = {
    "random_DM_phase": False,
    "random_DM_amplitude": False,
    "stochasticity_of_DM_signal_neglected": True,
}
SUPPORTED_POWER_DOMAINS = {
    "NATIVE_PS_rad2_after_C_or_idealized_synthetic_PS_NOT_TIME_DOMAIN",
    "AFTER_NATIVE_C_NOISE_PS_rad2_NOT_TIME_DOMAIN_DM",
}

def validate_paper_signal_simulation_manifest(manifest):
    if not isinstance(manifest, dict):
        raise ValueError("Simulation manifest must be a dictionary")
    for field, value in REQUIRED_FIXED_FLAGS.items():
        if manifest.get(field) is not value:
            raise ValueError(f"2026 paper approximation requires {field}={value}")
    if manifest.get("random_DM_velocity_modes", False) is not False:
        raise ValueError("Never sample random dark-matter velocity modes for deterministic template")
    if manifest.get("injection_domain") not in SUPPORTED_POWER_DOMAINS:
        raise ValueError("A deterministic axion POWER injection domain must be explicit")
    if manifest.get("source_paper_doi") != "10.1038/s41598-025-33428-2":
        raise ValueError("Source paper provenance mismatch")
    return True
