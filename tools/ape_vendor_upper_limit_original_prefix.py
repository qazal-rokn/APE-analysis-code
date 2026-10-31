"""Direct original vendor upper-limit numerical function, no AST extraction.

The original calc_upper_limit.main cannot be called unmodified for APE:
it unconditionally reads LIGO files and L1 transfer functions. Its portable
numerical calculation is get_upper_limits_approx, called here unchanged.
"""
from __future__ import annotations
import numpy as np
from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
from scalar_dark_matter_lpsd.lvk.calc_upper_limit import get_upper_limits_approx


def run_vendor_upper_limit_prefix(max_lkl_path, output_file, pruning=1,
                                  alpha_CL=.95, tol=1e-10, kernel_size=100000):
    max_lkl = MaxLKLIO(infile=str(max_lkl_path))
    if pruning < 1:
        raise ValueError("pruning must be >= 1")
    result = np.zeros((len(max_lkl.frequency[::pruning]), 4))
    result[:, 0] = max_lkl.frequency[::pruning]
    for i, sigma in enumerate(max_lkl.sigma[::pruning]):
        result[i, 1:3] = get_upper_limits_approx(
            sigma, alpha_CL=alpha_CL, tol=tol)
    np.savetxt(output_file, result, fmt="%.6e", delimiter="\t",
               header="Frequency (Hz)\tUpper limit (GeV^-1)\tUncertainty (GeV^-1)\tLIGO-only conversion omitted")
    return result
