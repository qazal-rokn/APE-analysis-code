# Verified segment extraction
Uses only the approved 0011/0012/0013 slices. Each source is processed from raw
index zero through its declared stop with unchanged reference quality/calibration
mapping; only then is its declared slice retained. Both recorded quadratures must
report unit `V`; every retained sample must be sample-good and calibration-valid.
No cache, HDF5 export, C execution, or spectrum is involved.
