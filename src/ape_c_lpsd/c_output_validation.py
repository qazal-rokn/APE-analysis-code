from __future__ import annotations
import math,re
import numpy as np
class COutputValidationError(RuntimeError): pass
def frequency_tolerance(token:str, parsed:float)->float:
 m=re.fullmatch(r'[+-]?(?:\d+)(?:\.(\d+))?[eE]([+-]?\d+)',token)
 if not m: raise COutputValidationError(f'frequency token is not scientific notation: {token!r}')
 decimals=len(m.group(1) or ''); exponent=int(m.group(2)); rounding=0.5*10.0**(exponent-decimals); allowance=4*np.spacing(abs(parsed))
 return rounding+allowance
def validate_frequency_rows(lines,expected):
 rows=[]; tol=[]; dif=[]
 for line in lines:
  if line.strip() and not line.startswith('#'):
   v=line.split()
   if len(v)!=5: raise COutputValidationError('expected exactly five columns')
   rows.append([float(x) for x in v]);tol.append(frequency_tolerance(v[0],float(v[0])))
 a=np.asarray(rows,float); expected=np.asarray(expected,float);tol=np.asarray(tol)
 if a.shape!=(len(expected),5) or not np.all(np.isfinite(a)) or not np.all(np.diff(a[:,0])>0) or not np.all(a[:,1]>=0) or not np.all(a[:,2]>0): raise COutputValidationError('invalid row count/order/PSD/average values')
 dif=np.abs(a[:,0]-expected)
 if np.any(dif>tol):
  i=int(np.flatnonzero(dif>tol)[0]);raise COutputValidationError(f'frequency mismatch row {i}: printed={a[i,0]} expected={expected[i]} difference={dif[i]} tolerance={tol[i]}')
 return a,{'maximum_observed_rounding_difference_hz':float(dif.max()),'per_row_frequency_tolerance_hz':tol.tolist()}
