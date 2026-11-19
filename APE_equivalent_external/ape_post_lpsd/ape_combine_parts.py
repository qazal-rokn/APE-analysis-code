"""Small APE provenance guard around the UNCHANGED vendor MLE part combiner.

The original combine_maxlkl_output.main() handles all numerical combining.
This APE adapter only verifies X/Y, response, and contiguous chunk identity.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path

import numpy as np

from scalar_dark_matter_lpsd.lvk.combine_maxlkl_output import main as vendor_combine
from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO


class APECombineError(ValueError):
    pass


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def combine_ape_parts(results_dir: str | Path, prefix: str, output_npy: str | Path):
    directory = Path(results_dir).expanduser().resolve()
    output = Path(output_npy).expanduser().resolve()
    provenance_path = output.with_suffix('.json')
    if (not directory.is_dir() or not isinstance(prefix, str) or not prefix
            or '/' in prefix or '\\' in prefix or any(c in prefix for c in '*?[]')):
        raise APECombineError('Invalid APE MLE part directory or prefix')
    if output.suffix != '.npy' or not output.parent.is_dir() or output.exists() or provenance_path.exists():
        raise APECombineError('Expected new .npy output and unused .json provenance in an existing folder')
    parts = sorted(Path(x) for x in glob.glob(str(directory / f'{prefix}_*_result.npy')))
    if not parts:
        raise APECombineError('No original-format Step-4 likelihood part outputs')
    manifest_list = []
    for part in parts:
        sibling = part.with_suffix('.json')
        if not sibling.is_file():
            raise APECombineError(f'Missing Step-4 manifest for {part.name}')
        with sibling.open(encoding='utf-8') as f:
            manifest = json.load(f)
        if manifest.get('channel') not in ('X','Y','Q') or manifest.get('joint_xy_likelihood_used') is not False:
            raise APECombineError('Only separate APE X/Y MLE parts may be combined')
        indices = manifest.get('native_grid_indices')
        if (not isinstance(indices, list) or len(indices) != 2 or
                any(not isinstance(i, int) for i in indices) or indices[0] >= indices[1]):
            raise APECombineError('Step-4 MLE parts must declare native grid chunk indices')
        values = MaxLKLIO(infile=str(part)).data
        if len(values) == 0 or values.shape[1] != 5 or not np.isfinite(values[:,0]).all():
            raise APECombineError('Invalid original five-column MLE part')
        if not np.all(np.diff(values[:,0]) > 0):
            raise APECombineError('MLE parts must already be internally frequency ordered')
        if manifest.get('frequency_count') != len(values):
            raise APECombineError('MLE part manifest row count mismatch')
        manifest_list.append((indices, part, sibling, manifest))
    manifest_list.sort(key=lambda x: x[0][0])
    first = manifest_list[0][3]
    base = manifest_list[0][0][0] - first.get('start_index_within_analysis_range', -1)
    if first.get('start_index_within_analysis_range') != 0 or base < 0:
        raise APECombineError('A full APE run must begin at the first requested analysis bin')
    for i, (indices, part, sibling, manifest) in enumerate(manifest_list):
        for key in ('channel', 'source_spectrum','source_background','source_response',
                    'response_metadata','signal_model','frequency_range_hz','scan_log10mu'):
            if manifest.get(key) != first.get(key):
                raise APECombineError(f'Inconsistent APE MLE chunk: {key}')
        if manifest.get('joint_xy_likelihood_used') is not False:
            raise APECombineError('Invalid X/Y statistical assumption')
        if indices[0] - manifest.get('start_index_within_analysis_range', -1) != base:
            raise APECombineError('Chunk start offset disagrees with native-grid indices')
        if i and indices[0] != manifest_list[i-1][0][1]:
            raise APECombineError('MLE chunk indices have a gap or overlap')
    try:
        # The original vendor is the only code that concatenates and sorts
        # the five-column MLE numerical results.
        vendor_combine(results_dir=str(directory), prefix=prefix, output_file=str(output), verbose=False)
        combined = MaxLKLIO(infile=str(output)).data
        expected_rows = sum(x[3]['frequency_count'] for x in manifest_list)
        if (len(combined) != expected_rows or
                np.any(np.diff(combined[:,0]) <= 0) or not np.isfinite(combined[:,0]).all()):
            raise APECombineError('Original combined output contains duplicate/invalid frequencies')
        report = dict(first)
        report['combined_parts'] = [
            {'part': str(part), 'sha256': _hash(part), 'manifest': str(sibling),
             'manifest_sha256': _hash(sibling), 'native_grid_indices': indices}
            for indices, part, sibling, _ in manifest_list]
        report['frequency_count'] = len(combined)
        report['native_grid_indices'] = [manifest_list[0][0][0], manifest_list[-1][0][1]]
        report['status'] = 'CONDITIONAL: original MLE parts combined; APE physics and coverage not certified'
        with provenance_path.open('x', encoding='utf-8') as f:
            json.dump(report, f, indent=2)
            f.write('\n')
    except BaseException:
        output.unlink(missing_ok=True)
        provenance_path.unlink(missing_ok=True)
        raise
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir', required=True)
    p.add_argument('--prefix', required=True)
    p.add_argument('--output-npy', required=True)
    a = p.parse_args()
    report = combine_ape_parts(a.results_dir, a.prefix, a.output_npy)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
