"""Adapter that calls, rather than reimplements, frozen APE Modules 2A, 2B, and 3."""
from __future__ import annotations
from pathlib import Path
import sys
from typing import Any
from .reference_paths import resolve_reference_dewesoft_paths

class ReferenceWorkflowError(RuntimeError): pass

REFERENCE_COMPONENTS=(
 "ape_analysis.quality.cavity_lock.build_cavity_lock_mask",
 "ape_analysis.quality.cavity_lock.CavityLockParameters",
 "ape_analysis.quality.pem_modulation.build_pem_modulation_quality",
 "ape_analysis.quality.pem_modulation.PemModulationParameters",
 "ape_analysis.calibration.polarization.build_polarization_calibration",
 "ape_analysis.calibration.polarization.PolarizationCalibrationParameters",
)

class ReferenceWorkflow:
    """Load frozen reference code from a configured read-only project path."""
    def __init__(self, reference_project: Path, cleaning_settings: dict[str,Any]):
        self.dewesoft_paths=resolve_reference_dewesoft_paths(reference_project)
        source=reference_project/"src"
        if not source.is_dir(): raise ReferenceWorkflowError(f"reference src missing: {source}")
        if str(source) not in sys.path: sys.path.insert(0,str(source))
        try:
            from ape_analysis.quality import CavityLockParameters, ExternallyReducedCavityRecord, PemModulationParameters, ReducedAverageChannel, build_cavity_lock_mask, build_pem_modulation_quality
            from ape_analysis.calibration import PolarizationCalibrationParameters, RawPemQuadratureChannel, build_polarization_calibration
            from ape_analysis.io import DewesoftReducedDatasetSource
        except ImportError as e: raise ReferenceWorkflowError("unable to import configured frozen reference workflow") from e
        self.CavityLockParameters=CavityLockParameters; self.ExternallyReducedCavityRecord=ExternallyReducedCavityRecord; self.PemModulationParameters=PemModulationParameters; self.ReducedAverageChannel=ReducedAverageChannel; self.build_cavity_lock_mask=build_cavity_lock_mask; self.build_pem_modulation_quality=build_pem_modulation_quality; self.PolarizationCalibrationParameters=PolarizationCalibrationParameters; self.RawPemQuadratureChannel=RawPemQuadratureChannel; self.build_polarization_calibration=build_polarization_calibration; self.DewesoftReducedDatasetSource=DewesoftReducedDatasetSource
        self.cleaning_settings=cleaning_settings; self._validate_settings_match_reference_defaults()

    def _validate_settings_match_reference_defaults(self) -> None:
        expected={"cavity_lock":{"f_lock":.90,"reference_window_seconds":60.0,"reference_quantile":.90,"refinement_iterations":3},"pem_modulation":{"cavity_transmission_mw_per_v":20.0,"second_harmonic_peak_filter_factor":.7,"responsivity_a_per_w":.7,"transimpedance_gain_v_per_a":1.e6,"rms_to_peak_factor":2**.5,"eta0_min_accepted":.019},"polarization":{"sample_rate_hz":20000.0,"reduced_block_size_seconds":.05,"samples_per_reduced_block":1000,"responsivity_a_per_w":.7,"transimpedance_gain_v_per_a":1.e6,"first_harmonic_rms_factor":2**.5}}
        for group,defaults in expected.items():
            actual=self.cleaning_settings.get(group)
            if actual != defaults: raise ReferenceWorkflowError(f"{group} differs from frozen reference defaults; approval/provenance review required")

    def quality_for_file_prefix(self, path: Path, reduced_count: int):
        """Exact reference quality functions; no cleaning rule is implemented here."""
        def read(name): return self.DewesoftReducedDatasetSource(path,reduced_record_start=0,reduced_record_count=reduced_count,**self.dewesoft_paths.source_kwargs()).read_reduced_channel(name)
        cavity,x,y=(read(n) for n in ("cavity_trans","pem2_X","pem2_Y")); origin=str(path); start=float(cavity.metadata["start_measure_time"])
        lock=self.build_cavity_lock_mask(tuple(self.ExternallyReducedCavityRecord(timestamp=r.dewesoft_timestamp,cavity_average=r.ave,cavity_minimum=r.min,source_file=origin,start_measure_time=start) for r in cavity.records),self.CavityLockParameters())
        def channel(name,v): return self.ReducedAverageChannel(recorded_name=name,timestamps=tuple(r.dewesoft_timestamp for r in v.records),averages=tuple(r.ave for r in v.records),block_size_seconds=v.block_size_seconds,source_file=origin,start_measure_time=start)
        return self.build_pem_modulation_quality(channel("pem2_X",x),channel("pem2_Y",y),channel("cavity_trans",cavity),lock.cavity_good,self.PemModulationParameters())
