"""Configuration A: call the ORIGINAL single-bin vendor optimizer through
the existing minimal APE-only response adapter. No halo/Monte-Carlo imports.
"""
from __future__ import annotations

from ape_post_lpsd.ape_max_likelihood import main as _vendor_adapter_main


if __name__ == "__main__":
    _vendor_adapter_main()
