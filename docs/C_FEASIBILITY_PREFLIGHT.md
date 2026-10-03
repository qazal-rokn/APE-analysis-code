# Strict C-LPSD Feasibility Preflight

Before HDF5 export or C execution, the bridge mirrors the C source's requested
node calculation: for each logarithmic node it computes
`fres=f*(exp(log(fmax/fmin)/(Jdes-1))-1)` and `required=round(fs/fres)`.
If `required` exceeds available uninterrupted samples, the run refuses before
data export or `lpsd-exec`. Its error identifies the first failing frequency,
node index, Jdes, required/available samples, and effective window settings.

The preflight never drops nodes, reduces density, or changes limits/window/PSLL/
overlap/averaging. A future real run must use it; if Q=1e6 is infeasible for the
globally longest segment, the researcher must make a new scientific decision.
