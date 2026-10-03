# Analysis Run Record Format

Status: proposed schema; populate only after contract approval.

Each run record must include a unique run ID; UTC timestamps; project revision;
contract/approval identifiers; input metadata and hashes; selected quantity and
unit; full segment inventory and selected endpoints; sample rate/count; HDF5
dataset name/dtype/shape/hash; executable and C-source hashes; build/platform and
dynamic-library details; exact controlled environment and argument vector; every
effective LPSD option; output/log/artifact hashes; validation IDs/results; exit
status/runtime; warnings; and the explicit real-data authorization when relevant.

Records are append-only artifacts. A rerun creates a new record.

