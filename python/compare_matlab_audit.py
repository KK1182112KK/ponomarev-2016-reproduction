"""Re-run the unchanged Python cases declared by the MATLAB audit artifact.
Usage: python python/compare_matlab_audit.py --matlab-dir /path/to/extracted/artifact
The script rejects mismatched source hashes and mismatched time grids. It does
not infer a successful theorem reproduction from cross-language agreement.
"""
from pathlib import Path
import argparse, importlib, hashlib, json, sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python/src'))
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--matlab-dir',type=Path,required=True)
ap.add_argument('--out',type=Path,default=root/'results/cross-language-audit')
a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
for line in (a.matlab_dir/'source-sha256.txt').read_text().splitlines():
    expected,name=line.split(None,1);p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
        raise RuntimeError('Source differs from MATLAB tested checkout: '+name)
meta=json.loads((a.matlab_dir/'audit.json').read_text())
functions={'scalar':('example_a','run_example_a'),'cascade':('example_b','run_example_b'),
 'pendulum':('example_c','run_example_c'),'predictor':('inexact_predictor','run_inexact_predictor'),
 'uncompensated':('uncompensated','run_uncompensated'),'delay-free':('delay_free','run_delay_free')}
rows=[]
for case in meta['cases']:
    if case['status']!='completed':
        rows.append(dict(name=case['name'],status='matlab_case_failed',error=case['error']));continue
    module,function=functions[case['kind']]
    r=getattr(importlib.import_module(module),function)(case['parameters'])
    x=r['x_hist'];x=x[:,None] if x.ndim==1 else x
    u=r['u_hist'][:,None] if 'u_hist' in r else np.column_stack([r['u1_hist'],r['u2_hist']])
    columns=[r['t'],x,u]
    for key in ('z_hist','V_hist','y_hist','P_hist'):
        if key in r:columns.append(r[key])
    py=np.column_stack(columns);ma=np.loadtxt(a.matlab_dir/case['csv'],delimiter=',',ndmin=2)
    if py.shape!=ma.shape or not np.allclose(py[:,0],ma[:,0],rtol=0,atol=1e-14):
        raise RuntimeError('Shapes or time grids differ: '+case['name'])
    if not np.isfinite(py).all():raise RuntimeError('Nonfinite Python result: '+case['name'])
    np.savetxt(a.out/case['csv'],py,delimiter=',')
    ix=np.atleast_1d(case['columns']['physical_state']).astype(int)-1
    iu=np.atleast_1d(case['columns']['control']).astype(int)-1
    rows.append(dict(name=case['name'],source_equations=case['source_equations'],
        parameters=case['parameters'],max_componentwise_X_difference=float(np.abs(py[:,ix]-ma[:,ix]).max()),
        max_componentwise_U_difference=float(np.abs(py[:,iu]-ma[:,iu]).max()),
        final_norm_python=float(np.linalg.norm(x[-1])),final_norm_matlab=case['final_norm']))
record=dict(matlab_tested_checkout=meta['tested_checkout'],matlab_run_id=meta['run_id'],
    matlab_version=meta['matlab_version'],python=sys.version,all_source_hashes_match=True,cases=rows)
(a.out/'summary.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
