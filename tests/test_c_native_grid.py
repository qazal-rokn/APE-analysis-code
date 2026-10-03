from ape_c_lpsd.c_grid import plan_c_native_grid, c_frequency, c_required_dft_samples

def test_small_plan_is_dense_and_feasible():
    p=plan_c_native_grid(fmin_hz=10.,fmax_hz=100.,segment_samples=10_000,sample_rate_hz=1_000.)
    assert p.required_dft_samples.max() <= 10_000
    assert c_required_dft_samples(10., p.relative_step, 1_000.) <= 10_000
    assert c_required_dft_samples(10., __import__('math').exp(__import__('math').log(10.)/p.jdes)-1., 1_000.) > 10_000
    assert p.frequency_hz[0] == c_frequency(0,10.,100.,p.jdes)
    assert p.frequency_hz[-1] == c_frequency(p.jdes-1,10.,100.,p.jdes)
