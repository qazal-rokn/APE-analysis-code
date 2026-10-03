import numpy as np
from ape_c_lpsd.c_output_validation import validate_frequency_rows,COutputValidationError
def test_accepts_c_10_decimal_scientific_rounding():
 a,p=validate_frequency_rows(['6.5000064784e+01 1e-15 1 1e-4 1'],np.array([65.0000647842543]));assert p['maximum_observed_rounding_difference_hz']<=p['per_row_frequency_tolerance_hz'][0]
def test_rejects_shifted_frequency():
 try: validate_frequency_rows(['6.5001000000e+01 1e-15 1 1e-4 1'],np.array([65.0000647842543]))
 except COutputValidationError: return
 raise AssertionError('shifted node accepted')
