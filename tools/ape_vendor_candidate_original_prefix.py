"""Run the EXACT detector-independent prefix of original vendor get_candidates.main.

No copied q0, trial-factor, cleaning or clustering implementation. Parse the
installed vendor source AST and execute its original statements unchanged up
to the 'Gather LPSD data' boundary. Refuse changed vendor source layouts.

The omitted suffix is the H1/L1-only data manager, detector consistency,
debugger, and plots. The vendor cluster() is called by the original prefix.
EXPERIMENTAL / UNTESTED.
"""
from __future__ import annotations
import argparse
import ast
import inspect
import json
from pathlib import Path

import numpy as np


def run_vendor_candidate_prefix(max_lkl_path, output_csv, sigma_limit=5.0):
    from scalar_dark_matter_lpsd.lvk import get_candidates as vendor
    source = inspect.getsource(vendor.main)
    tree = ast.parse(source)
    fn = tree.body[0]
    if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
        raise RuntimeError("Unexpected vendor candidate entry point")
    # Keep EXACT original source statements from the top of main through
    # candidate_idx/candidate_lambda. The next original statement creates
    # DataManager, which loads H1/L1 transfer functions.
    stop = next((i for i, stmt in enumerate(fn.body)
                 if isinstance(stmt, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id == "global_mngr"
                         for t in stmt.targets)), None)
    if stop is None:
        raise RuntimeError("Vendor candidate main layout changed; fail closed")
    prefix = fn.body[:stop]
    rendered = ast.unparse(ast.Module(body=prefix, type_ignores=[]))
    required = ("MaxLKLIO(", "zero_lkl", "max_lkl", "norm.ppf(",
                "candidate_mask", "cluster(", "candidate_idx", "candidate_lambda")
    if not all(s in rendered for s in required):
        raise RuntimeError("Vendor candidate prefix missing required original operations")
    # The original code is executed, not translated, copied, or reimplemented.
    code = compile(ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[])),
        inspect.getsourcefile(vendor.main), "exec")
    scope = dict(vars(vendor))
    scope.update(max_lkl_path=str(Path(max_lkl_path).resolve(strict=True)),
                 sigma_limit=float(sigma_limit))
    exec(code, scope, scope)
    frequencies = np.asarray(scope["lkl"].frequency)[~scope["mask"]][scope["candidate_idx"]]
    output = Path(output_csv).resolve()
    if output.exists():
        raise FileExistsError(output)
    np.savetxt(output, frequencies, fmt="%.12e")
    return {
        "status": "EXPERIMENTAL_UNTESTED_NOT_VALIDATED",
        "candidate_stage": "EXACT original vendor get_candidates.main prefix, AST-executed",
        "original_vendor_cluster_function": True,
        "removed": "original H1/L1 consistency checks and plotting suffix",
        "n_clusters": len(frequencies),
        "candidate_csv": str(output),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--max-lkl-path", required=True)
    p.add_argument("--output-csv", required=True)
    p.add_argument("--sigma-limit", type=float, default=5.0)
    a = p.parse_args()
    print(json.dumps(run_vendor_candidate_prefix(a.max_lkl_path,
            a.output_csv, a.sigma_limit), indent=2))
