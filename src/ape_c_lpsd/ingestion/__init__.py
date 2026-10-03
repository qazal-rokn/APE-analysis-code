"""APE ingestion, including canonical raw .dxd reading and legacy review tools."""

from .config import (
    APEChannelNames,
    CalibrationConfig,
    DXDIngestionConfig,
    IngestionConfig,
    load_dxd_ingestion_config,
    load_ingestion_config,
)
from .dxd_reader import DXDReaderError, read_ape_dxd
from .models import APEDXDData
from .segments import (
    ValidRun,
    choose_global_longest,
    freeze_manual_segment,
    write_frozen_segment_manifest,
)
from .reference_workflow import ReferenceWorkflow

__all__ = [
    "APEChannelNames",
    "APEDXDData",
    "CalibrationConfig",
    "DXDIngestionConfig",
    "DXDReaderError",
    "IngestionConfig",
    "ReferenceWorkflow",
    "ValidRun",
    "choose_global_longest",
    "freeze_manual_segment",
    "load_dxd_ingestion_config",
    "load_ingestion_config",
    "read_ape_dxd",
    "write_frozen_segment_manifest",
]
