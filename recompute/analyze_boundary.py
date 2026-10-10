"""Use fresh basin/spectrum estimates to test held-out finite-time thresholds."""
from pathlib import Path
import csv, math, itertools, json, hashlib
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'outputs/pairwise_paper_revision_v2/third_batch/B_time_boundary'
res=ROOT/'outputs/pairwise_paper_revision_v2/third_batch/residence_feasibility'
import h5py
rows=list(csv.DictReader((P/'thresholds_summary.csv').open()))
assert [float(r['T']) for r in rows]==[50,55,60,65,70,75]
assert all(r['status']=='bracketed' and float(r['width'])<=1e-6 and int(r['sampledRecrossLower'])==0 and int(r['sampledRecrossUpper'])==0 for r in rows)
with h5py.File(res/'residence_prediction.mat') as f:
    out=f['out'];pb=float(out['Pb'][()].item());lam=float(out['lambda'][()].item())
    epsb=float(out['basinUncertainty'][()].item());epsl=float(out['lambdaUncertainty'][()].item())
cal=rows[:2];tbar=52.5
logamp=sum(math.log(float(r['midpoint'])-pb)+lam*float(r['T']) for r in cal)/2
out=[]
for r in rows:
    t=float(r['T']);obs=float(r['midpoint']);delta=math.exp(logamp-lam*t);prediction=pb+delta
    bounds=[]
    for b,l,c1,c2 in itertools.product([pb-epsb,pb+epsb],[lam-epsl,lam+epsl],
        [float(cal[0]['P0Lower']),float(cal[0]['P0Upper'])],[float(cal[1]['P0Lower']),float(cal[1]['P0Upper'])]):
        bounds.append(b+math.sqrt((c1-b)*(c2-b))*math.exp(-l*(t-tbar)))
    low,high=min(bounds),max(bounds);rel=abs(obs-prediction)/(obs-pb)
    total=max(abs(prediction-low),abs(high-prediction))+float(r['width'])/2
    out.append(dict(T=t,role='calibration' if t<=55 else 'holdout',observedP0=obs,observedOffset=obs-pb,
        predictedP0=prediction,predictedOffset=delta,relativeOffsetError=rel,predictionLower=low,predictionUpper=high,
        totalErrorHalfWidth=total,relativeTotalHalfWidth=total/(obs-pb),
        empiricalGatePassed=rel<=.1 and total/(obs-pb)<=.1))
with (P/'boundary_prediction.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
report=dict(status='fresh-empirical-boundary-prediction',lambdaIndependent=lam,Pb=pb,epsilonPb=epsb,epsilonLambda=epsl,
    calibrationTimes=[50,55],holdoutTimes=[60,65,70,75],
    inputSHA256=hashlib.sha256((P/'time_windows.mat').read_bytes()).hexdigest(),
    maximumHeldoutRelativeOffsetError=max(r['relativeOffsetError'] for r in out if r['role']=='holdout'),
    allEmpiricalGatesPassed=all(r['empiricalGatePassed'] for r in out),
    uncertaintyScope='Empirical Nx2001/Nx4001 refinement allowance, not a rigorous continuum error bound.')
(P/'boundary_prediction_protocol.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
