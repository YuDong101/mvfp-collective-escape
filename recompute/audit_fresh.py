"""Numerical and provenance checks for fresh Figures 1--5 and supporting theory."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import h5py
from audit_sources import audit_sources
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
report={'run_id':'cp_fresh_recompute_2026_10_10','old_calculated_data_imported':False,'uploaded':False}
report['included_figures']=[1,2,3,4,5]
report['excluded_from_submission']='Figure 6 numerical data; supplied code and settings support new calculations.'
assert read(ROOT/'logs/local_validation_complete.json')['status']=='local-results-validated'
for name in ['core_a_complete.json','core_b_complete.json']:
    assert (ROOT/'logs'/name).is_file()
assert (ROOT/'logs/auxiliary_complete.json').is_file()
report['source_integrity']=audit_sources(ROOT)
report['fig123']=read(ROOT/'data/fig123/validation.json');assert report['fig123']['passed']
folder=ROOT/'data/fig45';source=read(folder/'source_subset.json')
assert source['source_sha256']==sha(folder/'combined_source.mat')
report['fig4']=read(folder/'fig4_checks.json');report['fig5']=read(folder/'fig5_checks.json')
for key in ['fig4','fig5']:
    assert report[key]['source_sha256']==source['source_sha256']
assert report['fig4']['max_self_consistency_residual']<1e-10
assert max(report['fig4']['geometry_residual_max'].values())<1e-10
assert not report['fig4']['connecting_absent_branches']
assert report['fig5']['prediction_sha256']==sha(folder/'residence_prediction.csv')
assert report['fig5']['max_component_sum_residual']<1e-12
assert report['fig5']['max_representative_energy_increase']<1e-6
assert report['fig5']['displayed_active_rate_samples']>100
theory=ROOT/'outputs/pairwise_paper_revision_v2'
spectra=[]
for nx in [201,1001,2001,4001]:
    path=theory/f'second_batch/local_spectrum/spectrum_Nx{nx}.mat'
    with h5py.File(path) as f:
        s=f['s'];rate=float(s['lambda'][()].item());residual=float(s['rightResidualL1'][()].item())
        assert rate>0 and residual<1e-7 and bool(s['positiveModeVerified'][()].item())
        if nx==201:assert int(s['numericalPositiveCount'][()].item())==1
        spectra.append({'Nx':nx,'lambda':rate,'right_eigenpair_residual':residual,'sha256':sha(path)})
report['spectra']=spectra
res=theory/'third_batch/residence_feasibility'
rows=list(csv.DictReader((res/'residence_prediction.csv').open()))
assert len(rows)==12 and all(r['usable']=='1' and r['eventStatus']=='exited' for r in rows)
held=[r for r in rows if r['role']=='holdout'];assert len(held)==6 and all(r['heldoutPass']=='1' for r in held)
report['residence']={'trajectories':12,'held_out':6,'maximum_heldout_relative_error':max(float(r['relativeError']) for r in held),'passed':True}
report['boundary']=read(theory/'third_batch/B_time_boundary/boundary_prediction_protocol.json')
assert report['boundary']['allEmpiricalGatesPassed']
stems=['fig1_frozen_landscape_limitations','fig2_dynamic_feedback_mechanism','fig3_finite_time_readout_transition','fig4_self_consistent_branch_regimes','fig5_free_energy_basin_selection']
figures={}
for stem in stems:
    for ext in ['pdf','png','svg']:
        path=ROOT/'figures'/f'{stem}.{ext}';assert path.is_file() and path.stat().st_size>1000
        figures[path.name]=sha(path)
report['candidate_figures_sha256']=figures
report['status']='complete-fresh-figures-1-to-5-audit-passed'
report['scope']='Figures 1--5 and supporting local theory, including Nx2001/Nx4001 refinement and held-out predictions. Figure 6 is excluded; its stopped partial calculations are not treated as complete.'
(ROOT/'provenance/final_numerical_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','source_integrity','residence','boundary']},indent=2))
