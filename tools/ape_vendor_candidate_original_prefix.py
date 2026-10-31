"""Execute unchanged original vendor candidate-selection statements.

Only omit original get_candidates.main's LIGO-specific DataManager,
two-IFO consistency checks, interactive debugger, and plots.
Fail closed if the vendor source structure changes.
"""
from __future__ import annotations
import ast
import inspect
from pathlib import Path
import numpy as np


def run_vendor_candidate_prefix(max_lkl_path, output_csv, sigma_limit=5.):
    from scalar_dark_matter_lpsd.lvk import get_candidates as vendor
    tree = ast.parse(inspect.getsource(vendor.main))
    fn = tree.body[0]
    if not isinstance(fn, ast.FunctionDef) or fn.name != "main":
        raise RuntimeError("Unexpected original vendor candidate entry point")
    boundary = next((i for i, stmt in enumerate(fn.body)
                     if isinstance(stmt, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == "global_mngr"
                             for t in stmt.targets)), None)
    if boundary is None:
        raise RuntimeError("Original candidate DataManager boundary not found")
    prefix = fn.body[:boundary]
    original = ast.unparse(ast.Module(body=prefix, type_ignores=[]))
    for required in ("MaxLKLIO(", "norm.ppf(", "cluster(", "candidate_idx"):
        if required not in original:
            raise RuntimeError("Original candidate-selection source changed")
    # Preserve the ORIGINAL vendor candidate_freqs and candidate_lambdas
    # assignments from immediately after the LIGO-only DataManager setup.
    suffix = fn.body[boundary:]
    assignments = []
    for stmt in suffix:
        if isinstance(stmt, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id in
            {"candidate_freqs", "candidate_lambdas"} for t in stmt.targets
        ):
            assignments.append(stmt)
    if len(assignments) != 2:
        raise RuntimeError("Original candidate-frequency assignments changed")
    fn.body = prefix + assignments + [ast.Return(
        value=ast.Name(id="candidate_freqs", ctx=ast.Load()))]
    namespace = vars(vendor).copy()
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])),
                 inspect.getsourcefile(vendor.main), "exec"), namespace, namespace)
    # Original signature, no guessed keyword/argument substitutions.
    frequencies = namespace["main"](
        str(Path(max_lkl_path).resolve(strict=True)), None, None, None,
        "", None, None, sigma_limit=float(sigma_limit))
    dest = Path(output_csv)
    if dest.exists():
        raise FileExistsError(dest)
    # Original get_candidates.main output serialization, unchanged.\n    with open(dest, "w") as f:\n        f.write("\\n".join([f"{x:.6e}" for x in frequencies]))
    return {"candidate_csv": str(dest), "n_clusters": int(len(frequencies)),
            "source": "unchanged original get_candidates.main statistical prefix"}
