"""Run original vendor upper-limit main's lookup-table calculation unchanged.

Only the original LIGO-only plotting/strain-conversion suffix is omitted.
Original source AST is inspected at runtime and fails closed on changes.
"""
import ast
import inspect
import numpy as np


def run_vendor_upper_limit_prefix(max_lkl_path, output_file, pruning=1,
                                  alpha_CL=.95, tol=1e-10, kernel_size=100000):
    from scalar_dark_matter_lpsd.lvk import calc_upper_limit as vendor
    tree = ast.parse(inspect.getsource(vendor.main))
    fn = tree.body[0]
    if not isinstance(fn, ast.FunctionDef):
        raise RuntimeError("Unexpected vendor UL entry point")
    stop = next((i for i, statement in enumerate(fn.body)
                 if isinstance(statement, ast.Expr)
                 and isinstance(statement.value, ast.Call)
                 and isinstance(statement.value.func, ast.Attribute)
                 and isinstance(statement.value.func.value, ast.Name)
                 and statement.value.func.value.id == "plt"
                 and statement.value.func.attr == "figure"), None)
    if stop is None:
        raise RuntimeError("Original vendor UL lookup-table boundary changed")
    prefix = fn.body[:stop]
    original = ast.unparse(ast.Module(body=prefix, type_ignores=[]))
    for operation in ("MaxLKLIO(", "np.geomspace(", "get_upper_limits_approx(",
                      "interp1d(", "ul_table(", "sigma_table("):
        if operation not in original:
            raise RuntimeError("Original vendor UL statistical workflow changed")
    fn.body = prefix + [ast.Return(value=ast.Name(id="result", ctx=ast.Load()))]
    namespace = dict(vars(vendor))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])),
                 inspect.getsourcefile(vendor.main), "exec"), namespace, namespace)
    result = namespace["main"](str(max_lkl_path), str(output_file),
                               pruning=pruning, alpha_CL=alpha_CL,
                               tol=tol, kernel_size=kernel_size)
    np.savetxt(output_file, result, delimiter="\t",
               header="Frequency (Hz)\tUpper limit (GeV^-1)\tUncertainty (GeV^-1)\tLIGO-only conversion omitted")
    return result
