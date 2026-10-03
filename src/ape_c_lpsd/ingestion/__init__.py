"""APE ingestion, including canonical raw .dxd reading and legacy review tools."""

from .config import (
    APEChannelNames,
    CalibrationConfig,
    DXDIngestionConfig,
    IngestionConfig,
    load_dxd_ingestion_config,
    load_ingestion_config,
)
from .dxd_reader import (
    DEFAULT_CHUNK_SIZE_SAMPLES,
    DXDReaderError,
    inspect_dxd,
    iter_ape_dxd_chunks,
)
from .models import APEDXDChunk, APEDXDData, APEDXDInspection
from .segments import (
    ValidRun,
    choose_global_longest,
    freeze_manual_segment,
    write_frozen_segment_manifest,
)
from .reference_workflow import ReferenceWorkflow

__all__ = [
    "APEChannelNames",
    "APEDXDChunk",
    "APEDXDData",
    "APEDXDInspection",
    "CalibrationConfig",
    "DEFAULT_CHUNK_SIZE_SAMPLES",
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
    "write_frozen_segment_manifest",
]
