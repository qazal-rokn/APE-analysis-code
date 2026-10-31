"""Execute original vendor UL lookup-table statements without LIGO-only suffix.

Original calc_upper_limit.main builds a 1000-point sigma lookup table,
interpolates, then performs hardcoded LIGO comparisons/plots and L1
conversion. Keep the original numerical prefix as AST, unchanged.
"""
from __future__ import annotations
import ast
import inspect
from pathlib import Path
import numpy as np


def run_vendor_upper_limit_prefix(max_lkl_path, output_file, pruning=1,
                                  alpha_CL=.95, tol=1e-10, kernel_size=100000):
    from scalar_dark_matter_lpsd.lvk import calc_upper_limit as vendor
    tree = ast.parse(inspect.getsource(vendor.main))
    fn = tree.body[0]
    if not isinstance(fn, ast.FunctionDef) or fn.name != "main":
        raise RuntimeError("Unexpected vendor upper-limit main")
    # Original statement '# Plot & save results' starts with plt.figure().
    boundary = next((i for i, stmt in enumerate(fn.body)
                     if isinstance(stmt, ast.Expr)
                     and isinstance(stmt.value, ast.Call)
                     and isinstance(stmt.value.func, ast.Attribute)
                     and isinstance(stmt.value.func.value, ast.Name)
                     and stmt.value.func.value.id == "plt"
                     and stmt.value.func.attr == "figure"), None)
    if boundary is None:
        raise RuntimeError("Original vendor upper-limit plotting boundary changed")
    prefix = fn.body[:boundary]
    original = ast.unparse(ast.Module(body=prefix, type_ignores=[]))
    for required in ("MaxLKLIO(", "np.geomspace(", "get_upper_limits_approx(",
                     "interp1d(", "ul_table(", "sigma_table("):
        if required not in original:
            raise RuntimeError("Original vendor upper-limit lookup-table code changed")
    fn.body = prefix + [ast.Return(value=ast.Name(id="result", ctx=ast.Load()))]
    namespace = vars(vendor).copy()
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])),
                 inspect.getsourcefile(vendor.main), "exec"), namespace, namespace)
    result = namespace["main"](
        str(Path(max_lkl_path).resolve(strict=True)), str(output_file),
        pruning=int(pruning), alpha_CL=float(alpha_CL), tol=float(tol),
        plot_path=None, kernel_size=float(kernel_size))
    # Original vendor result columns 1-2 are unchanged lookup-table outputs.
    # Column 3 is deliberately not filled: L1-strain equivalent is LIGO-only.
    dest = Path(output_file)
    if dest.exists():
        raise FileExistsError(dest)
    np.savetxt(dest, result[:, :3], fmt="%.6e", delimiter="\t",
               header="Frequency (Hz)\tUpper limit (GeV^-1)\tUncertainty (GeV^-1)")
    return result
