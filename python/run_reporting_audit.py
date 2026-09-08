"""Archive executions of the EXISTING solvers without changing their equations.
Run from any directory: python python/run_reporting_audit.py
Successful execution is not a stability proof or a paper-figure match.
"""
from __future__ import annotations
import contextlib, hashlib, io, json, platform, subprocess, sys, time, traceback
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python/src'))
from example_a import run_example_a
from example_b import run_example_b
from example_c import run_example_c

def main():
    out = ROOT / 'results/reporting-audit'
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    cases = [('scalar',run_example_a,{'dt':d}) for d in (.002,.001)]
    cases += [('cascade',run_example_b,{'dt':d}) for d in (.002,.001)]
    cases += [('pendulum',run_example_c,{'x0':[x,np.pi/2],'dt':.01}) for x in (np.pi/2,np.pi,3*np.pi/2)]
    cases += [('pendulum-refinement',run_example_c,{'dt':d}) for d in (.02,.005)]
    for i,(name,fn,params) in enumerate(cases):
        row={'case':name,'function':fn.__name__,'overrides':params,'status':'pending'}
        start=time.perf_counter(); log=io.StringIO()
        try:
            with contextlib.redirect_stdout(log):
                result=fn(params)
            t=np.asarray(result['t']); x=np.asarray(result['x_hist']).reshape(len(t),-1)
            u=np.asarray(result['u_hist']).reshape(len(t),-1)
            finite=bool(np.isfinite(np.column_stack((x,u))).all())
            row.update(status='completed' if finite else 'nonfinite',finite=finite,
                       t_end_actual=float(t[-1]),sample_count=len(t),initial_X=x[0].tolist(),
                       final_X=x[-1].tolist(),final_X_norm=float(np.linalg.norm(x[-1])),
                       initial_U=u[0].tolist(),max_abs_U=float(np.nanmax(np.abs(u))),
                       max_X_norm=float(np.nanmax(np.linalg.norm(x,axis=1))))
            path=f'{i:02d}-{name}.csv'
            np.savetxt(out/path,np.column_stack((t,x,u)),delimiter=',',comments='',
                       header=','.join(['t']+[f'X{j+1}' for j in range(x.shape[1])]+['U']))
            row['trajectory']=path
            if 'V_hist' in result: row['final_V']=float(result['V_hist'][-1])
        except Exception:
            row.update(status='exception',traceback=traceback.format_exc())
        row['elapsed_seconds']=time.perf_counter()-start
        (out/f'{i:02d}-{name}.log').write_text(log.getvalue(),encoding='utf-8')
        rows.append(row)
        (out/'runs.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(row),flush=True)
    def git(*args):
        p=subprocess.run(['git','-C',str(ROOT),*args],capture_output=True,text=True)
        return p.stdout.strip() if p.returncode==0 else None
    metadata={'created_utc':datetime.now(timezone.utc).isoformat(),'tested_checkout':git('rev-parse','HEAD'),
              'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted((ROOT/'python/src').glob('*.py'))},
              'scope':'Existing physical-plant solvers; no solver edits; all declared cases retained.',
              'limitations':['Nearest-grid history lookup and existing guards are preserved.',
                             'No pixel-based paper comparison or theorem certificate.']}
    (out/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    if any(r['status']!='completed' for r in rows): raise SystemExit(1)

if __name__=='__main__': main()
