"""APE ingestion, including canonical raw .dxd reading and legacy review tools."""

from .config import (
    APEChannelNames,
    CalibrationConfig,
    CleaningConfig,
    DXDIngestionConfig,
    IngestionConfig,
    load_dxd_ingestion_config,
    load_ingestion_config,
)
from .dxd_reader import (
    DXDReaderError,
    inspect_dxd,
    iter_ape_dxd_chunks,
    read_ape_dxd,
)
from .models import APECalibratedChunk, APEDXDChunk, APEDXDData, APEDXDInspection, APEDXDMetadata
from .segments import (
    ValidRun,
    choose_global_longest,
    freeze_manual_segment,
    write_frozen_segment_manifest,
)
from .reference_workflow import ReferenceWorkflow

__all__ = [
    "APECalibratedChunk",
    "APEChannelNames",
    "APEDXDChunk",
    "APEDXDData",
    "APEDXDInspection",
    "APEDXDMetadata",
    "CalibrationConfig",
    "CleaningConfig",
    "DXDIngestionConfig",
    "DXDReaderError",
    "IngestionConfig",
    "ReferenceWorkflow",
    "ValidRun",
    "choose_global_longest",
    "freeze_manual_segment",
    "inspect_dxd",
    "iter_ape_dxd_chunks",
    "load_dxd_ingestion_config",
    "load_ingestion_config",
    "read_ape_dxd",
    "write_frozen_segment_manifest",
]
