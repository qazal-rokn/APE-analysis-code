"""Explicit unboosted isotropic Maxwellian halo benchmark (not full SHM).

Uses nonrelativistic f = f_a*(1 + v^2/(2*c^2)). The user supplies
f_a, v0, c and bin edges. Lab velocity, escape cutoff and annual
modulation are NOT included; use a validated external halo distribution
for physical APE inference.
"""
import math
import numpy as np


def maxwell_halo_bin_weights(*, frequency_edges_hz, rest_frequency_hz,
                             velocity_scale_m_s, light_speed_m_s,
                             min_captured_probability):
    """Return frequency bin centers and normalized probability masses.

    Refuses under-covered grids rather than renormalizing away signal.
    """
    edges = np.asarray(frequency_edges_hz, dtype=float)
    p = [rest_frequency_hz, velocity_scale_m_s, light_speed_m_s,
         min_captured_probability]
    if not all(np.isscalar(x) and np.isfinite(x) for x in p):
        raise ValueError("Scalar halo inputs must be finite")
    if (rest_frequency_hz <= 0 or velocity_scale_m_s <= 0
            or light_speed_m_s <= 0 or velocity_scale_m_s >= light_speed_m_s
            or not 0 < min_captured_probability <= 1):
        raise ValueError("Invalid physical halo parameters")
    if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
            or np.any(np.diff(edges) <= 0) or edges[0] > rest_frequency_hz):
        raise ValueError("Edges must be increasing and cover rest frequency")
    v = light_speed_m_s*np.sqrt(
        2*np.maximum(edges/rest_frequency_hz - 1, 0))
    x = v/velocity_scale_m_s
    cdf = np.array([math.erf(t)-2/math.sqrt(math.pi)*t*math.exp(-t*t)
                    for t in x])
    weights = np.maximum(np.diff(cdf), 0)
    captured = weights.sum()
    if captured < min_captured_probability or captured <= 0:
        raise ValueError("Frequency grid misses too much halo probability")
    return (edges[:-1]+edges[1:])/2, weights/captured, captured
