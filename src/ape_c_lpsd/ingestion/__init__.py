"""Chronology-gated ingestion, review, and frozen-segment selection."""
from .config import IngestionConfig, load_ingestion_config
from .segments import ValidRun, choose_global_longest, freeze_manual_segment, write_frozen_segment_manifest
from .reference_workflow import ReferenceWorkflow

__all__ = ["IngestionConfig", "ReferenceWorkflow", "ValidRun", "choose_global_longest", "freeze_manual_segment", "load_ingestion_config", "write_frozen_segment_manifest"]
