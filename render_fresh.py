"""Render candidate figures from freshly computed data only."""
from pathlib import Path
import argparse, hashlib, importlib, json, sys

ROOT=Path(__file__).resolve().parent
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import MaxNLocator

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def figure5(module,out):
    source,folder=module._source(ROOT)
    dense=json.loads((folder/'dense_trajectories.json').read_text(encoding='utf-8'))
    selected=[next(c for c in dense['cases'] if c['P0']==p) for p in (0.,.071,.1)]
    critical=selected[1]
    source['freeEnergyBasinDiagnostics']['dissipationCase']=critical
    module._source=lambda root:(source,folder)
    module._representatives=lambda diagnostic:selected
    module.REMOVE_FITTED_FLUX_MARKER=False
    save=module.save_figure
    ranges={}
    def save_dense(fig,outdir,name):
        for ax in fig.axes[:2]:
            ax.set_xlim(0,100);ax.set_xticks([0,50,100])
            handles,labels=ax.get_legend_handles_labels()
            labels=['local' if label=='local potential' else label for label in labels]
            ax.legend(handles,labels,loc='center right',fontsize=8,handlelength=1.,labelspacing=.4)
        free=critical['freeEnergy'];tt=module.a(critical['t'])
        rate=module.a(free['minusDGdt']);rhs=module.a(free['rhsFreeEnergyDecay'])
        valid=np.isfinite(rate)&np.isfinite(rhs)&(rate>0)&(rhs>0)&np.asarray(free['dissipationComparisonMask'],bool)
        peak=np.argmax(np.where(valid&(tt>1),rate,-np.inf));active=valid&(rate>=1e-3*rate[peak]);lo=hi=peak
        while lo>0 and active[lo-1]:lo-=1
        while hi<len(tt)-1 and active[hi+1]:hi+=1
        pad=.05*(tt[hi]-tt[lo]);ax=fig.axes[2];ax.set_xlim(tt[lo]-pad,tt[hi]+pad)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3,integer=True))
        ax.set_ylim(rate[lo:hi+1].min()/1.4,max(rate[lo:hi+1].max(),rhs[lo:hi+1].max())*1.4)
        fig.axes[7].set_xlim(-.2,.3);fig.axes[7].set_xticks([-.2,0,.2])
        ranges.update(rate_display_time_window=[float(tt[lo]-pad),float(tt[hi]+pad)],
                      displayed_active_rate_samples=int(hi-lo+1),rate_peak_time=float(tt[peak]))
        save(fig,outdir,name)
    module.save_figure=save_dense
    checks=module.build_figure5(ROOT,out)
    checks.update(ranges)
    (folder/'fig5_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    return checks

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--figures',default='1,2,3,4,5');args=parser.parse_args()
    if any(n not in {1,2,3,4,5} for n in map(int,args.figures.split(','))):
        parser.error('This entry renders Figures 1--5; see README for Figure 6 parameter adaptations.')
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    report={}
    for n in map(int,args.figures.split(',')):
        print(f'Rendering fresh candidate figure {n}',flush=True)
        module=importlib.import_module('plotting.'+('fig123' if n<=3 else 'fig45' if n<=5 else 'fig6'))
        report[str(n)]=figure5(module,out) if n==5 else getattr(module,f'build_figure{n}')(ROOT,out)
    report['output_sha256']={p.name:sha(p) for p in out.iterdir() if p.suffix in {'.pdf','.png','.svg'}}
    report['approval_status']='candidate; author approval required before manuscript replacement'
    (ROOT/'provenance'/('render_'+args.figures.replace(',','_')+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')

if __name__=='__main__':main()
