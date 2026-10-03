#!/usr/bin/env python3
"""Plan the approved no-data C-native grid; run with --help."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import yaml
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from ape_c_lpsd.c_grid import plan_c_native_grid, write_manifest

def main() -> None:
    parser = argparse.ArgumentParser(description='Create an immutable C-native grid manifest without reading signals or running C.')
    parser.add_argument('--config', default=ROOT / 'configs/run.yaml', type=Path)
    parser.add_argument('--reference-grid', default=ROOT / 'reference_grid/final_native_65Hz_9kHz_exact_ccore_plan.npz', type=Path)
    args = parser.parse_args()
    with args.config.open() as handle: cfg = yaml.safe_load(handle)
    bounds = cfg['c_native_grid']
    seg = cfg['validated_segment']
    l = cfg['lpsd']
    sbin = float(bounds['effective_sbin_bins'])
    absolute_c_min = sbin / (int(seg['sample_count']) / float(seg['sample_rate_hz'])) * (1.0 + (1.0 - float(l['overlap_percent']) / 100.0) * (int(l['minimum_averages']) - 1))
    if float(bounds['fmin_hz']) * (1.0 + 1e-6) < absolute_c_min: raise RuntimeError(f'configured search fmin is below C absolute window bound: {absolute_c_min}')
    plan = plan_c_native_grid(fmin_hz=float(bounds['fmin_hz']), fmax_hz=float(bounds['fmax_hz']), segment_samples=int(seg['sample_count']), sample_rate_hz=float(seg['sample_rate_hz']))
    options = {'-h': l['method'], '-w': l['window'], '-p': l['psll_db'], '-r': l['linear_regression'], '-l': l['overlap_percent'], '-a': l['desired_averages'], '-m': l['minimum_averages'], '-x': l['input_scale'], '-M': l['max_memory_power'], '-s': plan.fmin_hz, '-t': plan.fmax_hz, '-J': plan.jdes, '-n': plan.jdes, '-N': 0, 'effective_window_sbin_bins': sbin, 'c_absolute_window_lower_bound_hz': absolute_c_min}
    path = write_manifest(plan=plan, output_directory=Path(cfg['output_directory']) / 'c-native-grid-manifests', c_options=options, reference_grid=args.reference_grid)
    print(json.dumps({'manifest': str(path), 'first_frequency_hz': float(plan.frequency_hz[0]), 'last_frequency_hz': float(plan.frequency_hz[-1]), 'node_count': plan.jdes, 'relative_step': plan.relative_step, 'maximum_required_dft_samples': int(plan.required_dft_samples.max()), 'c_absolute_window_lower_bound_hz': absolute_c_min, 'all_nodes_fit': bool((plan.required_dft_samples <= plan.segment_samples).all()), 'no_signal_data_read': True, 'no_hdf5_written': True, 'no_c_execution': True}, indent=2))
if __name__ == '__main__': main()
