"""Bounded-memory, lossless Constant-Q batch validation and combination."""
from __future__ import annotations
import math
from pathlib import Path
from .constantq import ConstantQError, ENDPOINTS, FREQUENCIES


def _rows(path, expected_count, *, columns):
    count = 0
    last = None
    with Path(path).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != columns:
                raise ConstantQError(f"{path}:{line_number}: expected {columns} columns")
            try:
                v = tuple(float(x) for x in fields)
            except ValueError as exc:
                raise ConstantQError(f"{path}:{line_number}: nonnumeric output") from exc
            if not all(math.isfinite(x) for x in v):
                raise ConstantQError(f"{path}:{line_number}: nonfinite output")
            f, psd, n_avg, n_kernel = v[:4]
            if f <= 0 or psd < 0 or n_avg < 0 or n_kernel < 0:
                raise ConstantQError(f"{path}:{line_number}: invalid output values")
            if last is not None and f <= last:
                raise ConstantQError(f"{path}:{line_number}: duplicate/descending frequency")
            last = f
            count += 1
            if count > expected_count:
                raise ConstantQError(f"{path}: too many bins")
            yield v
    if count != expected_count:
        raise ConstantQError(f"{path}: missing bins, expected {expected_count}, got {count}")


def validate_and_combine(paths, output, *, endpoints=ENDPOINTS, columns=6):
    """Validate exact indexed batches and write first four original columns.

    Does not resample, average, scale, smooth, or modify C PSD estimates.
    Atomic publication: invalid inputs never replace a finished output.
    """
    if len(paths) != len(endpoints)-1 or endpoints[0] != 0:
        raise ConstantQError("Wrong number of index batches")
    if any(b <= a for a, b in zip(endpoints, endpoints[1:])):
        raise ConstantQError("Invalid batch boundaries")
    if columns != 6:
        raise ConstantQError("Modern method-2 fDNbri must contain six columns")
    target = Path(output)
    if target.exists() or target.with_suffix(target.suffix+".partial").exists():
        raise ConstantQError("Refusing to overwrite existing/partial combination")
    tmp = target.with_suffix(target.suffix+".partial")
    previous = None
    written = 0
    try:
        with tmp.open("x", encoding="utf-8") as dst:
            for i, p in enumerate(paths):
                for row in _rows(p, endpoints[i+1]-endpoints[i], columns=columns):
                    if previous is not None and row[0] <= previous:
                        raise ConstantQError("Duplicate or misordered frequency at batch boundary")
                    previous = row[0]
                    dst.write(" ".join(format(x, ".17g") for x in row[:4])+"\n")
                    written += 1
        if written != endpoints[-1]:
            raise ConstantQError("Combined spectrum has wrong number of bins")
        tmp.rename(target)
        return written
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def check_xy_frequency_alignment(x_file, y_file):
    """Streaming comparison: X/Y files must contain the exact same native grid."""
    with Path(x_file).open() as x, Path(y_file).open() as y:
        for index, (a, b) in enumerate(zip(x, y), 1):
            if not a or not b or a.split()[0] != b.split()[0]:
                raise ConstantQError(f"X/Y grid mismatch at row {index}")
        if x.readline() or y.readline():
            raise ConstantQError("Different X/Y row counts")
