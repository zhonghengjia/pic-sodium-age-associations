"""Read R-generated designs, screen complete/quasi separation using bounded LP."""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.optimize import linprog
import scipy


def screen(matrix, outcome):
    assert np.isfinite(matrix).all() and set(np.unique(outcome)).issubset({0, 1})
    scales = matrix.std(axis=0)
    scales[scales < 1e-12] = 1
    if np.linalg.matrix_rank(matrix / scales) < matrix.shape[1]:
        return {'status': 'RANK_DEFICIENT', 'pass': False}
    z = matrix / scales * (2 * outcome - 1)[:, None]
    n, p = z.shape
    signed = np.column_stack((z, -z))
    objective = -signed.mean(axis=0)
    inequality = np.vstack((-signed, np.ones((1, 2 * p))))
    bounds = np.r_[np.zeros(n), 1.0]
    res = linprog(objective, A_ub=inequality, b_ub=bounds, bounds=(0, None),
                  method='highs', options={'primal_feasibility_tolerance': 1e-9,
                                           'dual_feasibility_tolerance': 1e-9})
    if not res.success:
        return {'status': 'LP_FAILED', 'message': res.message, 'pass': False}
    score = max(0.0, float(-res.fun))
    margins = signed @ res.x
    return {'status': 'SEPARATION_DETECTED' if score > 1e-8 else 'NO_SEPARATION_DETECTED_AT_TOLERANCE',
            'objective': score, 'minimum_margin': float(margins.min()),
            'objective_tolerance': 1e-8, 'solver_feasibility_tolerance': 1e-9,
            'column_scales': scales.tolist(), 'solver_status': int(res.status),
            'l1_norm_split_coefficients': float(res.x.sum()),
            'max_constraint_violation': float(max(0, -margins.min(), res.x.sum() - 1)),
            'pass': score <= 1e-8 and float(margins.min()) >= -1e-7}


def main(folder):
    folder = Path(folder)
    metadata = json.loads((folder / 'design_metadata.json').read_text(encoding='utf-8-sig'))
    results = []
    for row in metadata['designs']:
        n, p = row['rows'], row['columns']
        values = np.fromfile(folder / row['path'], dtype='<f8')
        assert len(values) == n * (p + 1)
        result = screen(values[n:].reshape(n, p), values[:n])
        results.append({'imputation': row['imputation'], 'n': n, 'p': p, **result})
    out = {'method': 'bounded linear-programming complete/quasi separation screen',
           'scipy': scipy.__version__, 'model': metadata['model'],
           'all_pass': all(r['pass'] for r in results), 'results': results}
    (folder / 'separation_validation.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps({'model': out['model'], 'designs': len(results), 'all_pass': out['all_pass']}))
    return 0 if out['all_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1]))
