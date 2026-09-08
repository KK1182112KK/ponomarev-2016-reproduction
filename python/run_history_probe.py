"""Record, do not hide, negative-time and quadrature endpoint discrepancies.
Run from a repository checkout: python python/run_history_probe.py
This is a diagnostic recorder, not a passing-test or theorem certificate.
"""
from pathlib import Path
import importlib, json, sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python/src'))
is_pono=(root/'python/src/example_a.py').exists()
mods=['example_a','example_b','example_c'] if is_pono else ['inexact_predictor','uncompensated']
rows=[]
for name in mods:
    m=importlib.import_module(name)
    buffer=np.r_[np.ones(50),9.,9.]
    actual=float(m._lookup_u(-.004,.01,buffer,50))
    rows.append(dict(module=name,query=-.004,dt=.01,negative_history=1.,issued_at_zero=9.,
                     expected_strict_prehistory=1.,actual=actual,strict_prehistory_preserved=actual==1.))
record=dict(lookup_probes=rows,scope='Synthetic non-grid query; not evidence that every default trajectory is wrong.')
if is_pono:
    from example_b import run_example_b
    r=run_example_b(dict(h=1.,dt=.001,t_end=2.,x0=[1.,1.],u0=0.))
    t,x,u,z=r['t'],r['x_hist'],r['u_hist'],r['z_hist']
    # Exact integral of the actual ZOH history, evaluated AFTER the run.
    prefix=np.r_[0.,.001*np.cumsum(u[:-1])]
    integral=np.array([prefix[k]-prefix[max(0,k-1000)] for k in range(len(t))])
    reconstructed=x[:,0]-x[:,1]+integral
    record['eq86_diagnostic']=dict(dt=.001,T=2.,max_error=float(np.abs(reconstructed-z[:,0]).max()),
        stored_z1_at_t1=float(z[1000,0]),zoh_z1_at_t1=float(reconstructed[1000]),
        interpretation='Legacy trapezoid weights an unissued zero endpoint; comparison is not used by the controller.')
out=root/'results/history-probe';out.mkdir(parents=True,exist_ok=True)
(out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
