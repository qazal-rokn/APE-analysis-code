# Quick Start

Use the project-local environment:

```sh
.venv/bin/python tools/ape_c_lpsd.py --help
.venv/bin/python tools/ape_c_lpsd.py inspect --config configs/run.yaml
.venv/bin/python tools/ape_c_lpsd.py synthetic --config configs/run.yaml
.venv/bin/python tools/run_synthetic_validation.py
```

`inspect` reads only cache metadata, file paths, and masks. `synthetic` creates
only temporary synthetic data; it does not open the configured cache.

The future command is deliberately blocked unless both approvals are supplied:

```sh
.venv/bin/python tools/ape_c_lpsd.py run-real --config configs/run.yaml --approve-real-data
```

It currently still refuses after those flags: a separate post-validation approval
and completed real-run implementation are required. Edit a copy of
`configs/run.example.yaml`, not source code. Every tool supports `--help`.
