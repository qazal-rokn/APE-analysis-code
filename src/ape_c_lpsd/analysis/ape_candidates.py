"""APE candidate preselection: vendor-inspired statistics, no detector veto.

The original H1/L1 consistency test is not portable to X/Y because they
are quadratures of one apparatus, not independent observatories. Candidate
p-values/global significance require null simulations and correlation study.
"""
import numpy as np
from scipy.stats import norm


def select_candidates(frequencies_hz, log_likelihood_null, log_likelihood_max,
                      mu_hat, *, local_sigma_threshold,
                      n_effective_trials=None):
    """Local one-sided q0 candidates; optional *assumed* trial-factor correction.

    Output is a list of candidate dicts, retaining the largest q0 in each
    contiguous above-threshold run; terminal clusters are not dropped.
    No physical detection significance or independent-detector veto claimed.
    """
    f=np.asarray(frequencies_hz,dtype=float)
    l0=np.asarray(log_likelihood_null,dtype=float)
    lm=np.asarray(log_likelihood_max,dtype=float)
    mu=np.asarray(mu_hat,dtype=float)
    if (f.ndim!=1 or not f.size or any(a.shape!=f.shape for a in (l0,lm,mu))
            or not all(np.isfinite(a).all() for a in (f,l0,lm,mu))
            or np.any(f<=0) or np.any(np.diff(f)<=0)):
        raise ValueError("Candidate arrays must be finite aligned frequency vectors")
    if not np.isscalar(local_sigma_threshold) or not np.isfinite(local_sigma_threshold) or local_sigma_threshold<=0:
        raise ValueError("Positive finite local sigma threshold required")
    local_p=norm.sf(float(local_sigma_threshold))
    if n_effective_trials is None:
        threshold=float(local_sigma_threshold)
        correction="none: local ranking only"
    else:
        if not isinstance(n_effective_trials,int) or n_effective_trials<1:
            raise ValueError("n_effective_trials must be positive integer")
        threshold=float(norm.isf(local_p/n_effective_trials))
        correction="user-assumed independent effective trials; not coverage calibrated"
    q0=np.where(mu>0,np.maximum(0,2*(lm-l0)),0)
    flagged=(q0>threshold**2)
    results=[]
    for i in range(len(f)):
        if not flagged[i]:
            continue
        if i and flagged[i-1]:
            continue
        j=i+1
        while j<len(f) and flagged[j]:
            j+=1
        best=i+int(np.argmax(q0[i:j]))
        results.append(dict(frequency_hz=float(f[best]),q0=float(q0[best]),
                            local_sqrt_q0=float(np.sqrt(q0[best])),
                            mu_hat=float(mu[best]),cluster_start_index=int(i),
                            cluster_end_index_exclusive=int(j),
                            status="candidate_for_followup_not_detection"))
    return dict(candidates=results,threshold_sigma=threshold,
                trial_correction=correction,total_bins=int(len(f)),
                flagged_bins=int(flagged.sum()))
