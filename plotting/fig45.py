"""Approved candidate redraw of branch geometry and free-energy evidence.

Figure 4: one root set determines branch-specific potentials and extrema.
Figure 5: fixed-functional trajectories select a basin; the finer local residence
prediction is explicitly a separate discretization. No dynamical solves here.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import h5py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import MaxNLocator, NullLocator

from .style import setup_style, panel_label, clean_axis, save_figure, COLORS


# Candidate-only export hook; the approved/default figure retains the fitted-flux marker.
REMOVE_FITTED_FLUX_MARKER = False


def a(value):
    return np.asarray(value, dtype=float).reshape(-1)


def _source(root):
    folder = Path(__file__).resolve().parents[1]/'data/fig45'
    return json.loads((folder/'source_subset.json').read_text(encoding='utf-8')), folder


def _line(color='black', ls='-', marker=None, label=''):
    return Line2D([], [], color=COLORS.get(color, color), ls=ls, marker=marker,
                  lw=1.3, ms=3.5, label=label)


def _manifest(root, folder, out):
    pass  # The candidate runner records all input and output hashes.


def _finish(ax, xlabel, ylabel):
    clean_axis(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))


def _branch_arrays(xb):
    """One position per scan value, with explicit NaNs at absent roots.

    Stable left/right mean negative/nonnegative M1 as in the source. Verify
    at most one root of each identity at each scan value before connecting.
    """
    sigma = a(xb['branchValue'])
    grid = np.unique(sigma)
    stable = np.asarray(xb['stable'], bool)
    mean = a(xb['xBranch'])
    masks = [stable & (mean < 0), stable & (mean >= 0), ~stable]
    output = []
    for mask in masks:
        ids = np.full(len(grid), -1, dtype=int)
        for j, value in enumerate(grid):
            same = np.flatnonzero(mask & (sigma == value))
            assert len(same) <= 1, 'Root identity requires continuation, not rank joining'
            if len(same):
                ids[j] = same[0]
        fields = {}
        for key in ('xBranch', 'F', 'Pss', 'xL', 'xs', 'xR', 'xSingle'):
            values = np.full(len(grid), np.nan)
            good = ids >= 0
            values[good] = a(xb[key])[ids[good]]
            fields[key] = values
        output.append(fields)
    return grid, output


def build_figure4(root: Path, out: Path):
    root, out = Path(root), Path(out)
    d, folder = _source(root)
    xb, par = d['xBifurcationVsSigma'], d['par']
    grid, groups = _branch_arrays(xb)
    setup_style()
    fig = plt.figure(figsize=(179/25.4, 130/25.4))
    gs = fig.add_gridspec(2, 3, left=.095, right=.975, top=.94, bottom=.12,
                         wspace=.48, hspace=.50, height_ratios=[1, 1.05])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    axes += [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1:])]
    colors = ['blue', 'red', 'purple']
    styles = ['-', '-', '--']
    for ax, field, ylabel in zip(axes[:3], ['xBranch', 'F', 'Pss'],
                               [r'$M_1$', r'$F^*$', r'$P_{\mathrm{ss}}$']):
        for group, color, ls in zip(groups, colors, styles):
            ax.plot(grid, group[field], color=COLORS[color], ls=ls)
        for value in [.5, 2]:
            ax.axvline(value, color='#aaaaaa', ls=':', lw=.9, zorder=0)
        ax.set_xscale('log')
        ax.set_xlim(grid[0], grid[-1])
        ax.set_xticks([.01, .1, 1, 10])
        _finish(ax, r'$\sigma_1$', ylabel)
        ax.xaxis.set_minor_locator(NullLocator())
    axes[2].set_ylim(-.035, 1.035)
    axes[2].set_yticks([0, .5, 1])
    handles = [_line(c, ls, label=label) for c, ls, label in zip(colors, styles,
               [r'stable left $F^*$', r'stable right $F^*$', r'unstable $F^*$'])]
    axes[1].legend(handles=handles, loc='upper left', handlelength=1.5,
                   fontsize=8.5,labelspacing=.35)

    sigma = a(xb['branchValue'])
    ix = np.flatnonzero(np.isclose(sigma, 2, atol=1e-12, rtol=0))
    ix = ix[np.argsort(a(xb['F'])[ix])]
    assert len(ix) == 3
    roots = a(xb['F'])[ix]
    assert list(np.asarray(xb['stable'])[ix]) == [True, False, True]
    xx = np.linspace(-1, 1, 1001)
    potential_rows = []
    for ff, color, ls in zip(roots, ['blue', 'purple', 'red'], ['-', '--', '-']):
        potential = -.5*(par['u']-2)*xx**2 + .25*(par['v']+par['sigma2'])*xx**4-ff*xx
        shifted = potential-potential.min()
        axes[3].plot(xx, shifted, color=COLORS[color], ls=ls,
                     label=fr'$F^*={ff:.3g}$')
        potential_rows.extend(zip(np.repeat(ff, len(xx)), xx, potential, shifted))
    axes[3].set_xlim(-1, 1)
    # Original view is a local potential window; values outside it are explicitly
    # retained in the exported source table and described in the caption.
    axes[3].set_ylim(-.03, 1)
    axes[3].text(.63,.50,r'$\sigma_1=2$',transform=axes[3].transAxes,
                 ha='center',va='center',fontsize=9)
    _finish(axes[3], r'$x$', 'Shifted potential')

    ax = axes[4]
    extrema_specs = [('xL', 'o', '-', 'left minimum'), ('xs', 'x', '--', 'saddle'),
                    ('xR', 's', '-', 'right minimum'), ('xSingle', 'D', ':', 'single minimum')]
    for group, color in zip(groups, colors):
        for key, marker, ls, label in extrema_specs:
            if np.any(np.isfinite(group[key])):
                ax.plot(grid, group[key], ls=ls, color=COLORS[color], marker=marker,
                        markevery=9, ms=3.5, lw=1.0, mfc='white' if marker!='x' else COLORS[color])
    ax.set_xscale('log')
    ax.set_xticks([.01, .1, 1, 10])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(grid[0], grid[-1])
    ax.set_ylim(-1, 1)
    ax.axvline(par['u'], color=COLORS['black'], ls='-.', lw=.9)
    ax.text(par['u']*.91, -.94, 'zero quadratic term', fontsize=8.5, ha='right', va='bottom')
    _finish(ax, r'$\sigma_1$', r'potential extrema $x$')
    extrema_handles=[_line('black',ls,marker,label) for _,marker,ls,label in extrema_specs]
    for handle in extrema_handles:
        if handle.get_marker()!='x': handle.set_markerfacecolor('white')
    ax.legend(handles=extrema_handles,
              loc='upper left', ncol=2, handlelength=1.5, columnspacing=1, labelspacing=.5)
    for label, ax in zip('abcde', axes):
        panel_label(ax, label)

    checks = {'source_sha256': d['source_sha256'], 'sigma_range': [grid[0], grid[-1]],
              'scan_values': len(grid), 'roots': len(sigma), 'F_at_exact_sigma1_2': roots.tolist(),
              'max_self_consistency_residual': float(np.nanmax(np.abs(a(xb['residualAtRoot'])))),
              'geometry_residual_max': {}, 'connecting_absent_branches': False,
              'potential_display_window': [-.03, 1], 'potential_data_x_range': [-1, 1]}
    for key in ('xL', 'xs', 'xR', 'xSingle'):
        x = a(xb[key]); valid = np.isfinite(x)
        residual = (par['v']+par['sigma2'])*x[valid]**3-(par['u']-sigma[valid])*x[valid]-a(xb['F'])[valid]
        checks['geometry_residual_max'][key] = float(np.max(np.abs(residual)))
    assert max(checks['geometry_residual_max'].values()) < 1e-8
    with (folder/'fig4_conditioned_potentials.csv').open('w', newline='', encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['Fstar','x','U','UminusMin']); w.writerows(potential_rows)
    (folder/'fig4_checks.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    save_figure(fig, out, 'fig4_self_consistent_branch_regimes')
    _manifest(root, folder, out)
    return checks


def _representatives(diag):
    cases = diag['cases']
    mid = float(a(diag['asymptoticBoundary']['midpoints'])[0])
    p = np.array([c['P0'] for c in cases])
    indices = [int(np.argmin(p)), int(np.argmax(p)),
               int(np.flatnonzero(p <= mid)[np.argmax(p[p <= mid])]),
               int(np.flatnonzero(p >= mid)[np.argmin(p[p >= mid])])]
    return [cases[i] for i in indices]


def build_figure5(root: Path, out: Path):
    root, out = Path(root), Path(out)
    d, folder = _source(root)
    diag = d['freeEnergyBasinDiagnostics']
    cases = sorted(diag['cases'], key=lambda c:c['P0'])
    reps = _representatives(diag)
    critical = diag['dissipationCase']
    predfile = folder/'residence_prediction.csv'
    with predfile.open(newline='', encoding='utf-8') as f:
        predictions = list(csv.DictReader(f))
    assert all(r['usable']=='1' and r['eventStatus']=='exited' for r in predictions)
    with h5py.File(predfile.with_suffix('.mat'),'r') as stored:
        local=stored['out']
        local_meta={k:float(local[k][()].item()) for k in
                    ('Pb','lambda','lambdaUncertainty','basinUncertainty')}
        local_meta.update({k:float(local['config'][k][()].item()) for k in
                           ('PbGridN','entryRadius','baseExitRadius','traceStep','T')})
    for row in predictions:
        assert abs(abs(float(row['P0'])-local_meta['Pb'])-float(row['distance'])) < 1e-14
    setup_style()
    fig = plt.figure(figsize=(179/25.4, 170/25.4))
    top = fig.add_gridspec(1, 3, left=.095, right=.955, bottom=.73, top=.94, wspace=.51)
    lower = fig.add_gridspec(2, 3, left=.095, right=.915, bottom=.09, top=.60,
                             wspace=.65,hspace=.48,width_ratios=[1, 1, 1.2],height_ratios=[1,1])
    axes = [fig.add_subplot(top[0,i]) for i in range(3)]
    axes += [fig.add_subplot(lower[0,i]) for i in range(2)]
    ftop=fig.add_subplot(lower[0,2]); fbottom=fig.add_subplot(lower[1,2]); axes.append(ftop)
    colors = ['purple', 'green', 'blue', 'orange']
    styles = ['-', '--', '-.', ':']
    for case, color, ls in zip(reps, colors, styles):
        axes[0].plot(a(case['t']), a(case['freeEnergy']['totalChange']),
                     color=COLORS[color], ls=ls, label=fr'$P_0={case["P0"]:.4g}$')
    axes[0].axhline(0, color='#999999', lw=.7, ls=':')
    _finish(axes[0], r'$t$', 'Free-energy change')
    handles=axes[0].get_legend_handles_labels()[0]
    axes[0].legend(handles=handles,loc='center left',bbox_to_anchor=(.03,.52),
                   handlelength=1.5,labelspacing=.4,fontsize=8.5)
    free=critical['freeEnergy']; tt=a(critical['t'])
    specs=[('total','black','-',r'$\Delta\mathcal{G}$'),('entropy','red','--','entropy'),
           ('localPotential','orange','-.','local potential'),('interaction','purple',':','interaction')]
    for key,color,ls,label in specs:
        yy=a(free[key]); axes[1].plot(tt,yy-yy[0],color=COLORS[color],ls=ls,label=label)
    _finish(axes[1],r'$t$','component change')
    axes[1].text(.05,.13,fr'$P_0={critical["P0"]:g}$',transform=axes[1].transAxes,
                 ha='left',fontsize=9)
    axes[1].legend(loc='center left',bbox_to_anchor=(.01,.48),ncol=1,labelspacing=.4,handlelength=1.5)
    minus,rhs=a(free['minusDGdt']),a(free['rhsFreeEnergyDecay'])
    valid=np.isfinite(minus)&np.isfinite(rhs)&(minus>0)&(rhs>0)&np.asarray(free['dissipationComparisonMask'],bool)
    active=valid&(minus>=1e-3*minus[valid].max())
    axes[2].plot(tt,np.where(active,minus,np.nan),color=COLORS['black'],label='Time diff.')
    axes[2].plot(tt,np.where(active,rhs,np.nan),color=COLORS['red'],ls='--',label='RHS rate')
    _finish(axes[2],r'$t$','energy decay rate')
    axes[2].set_yscale('log')
    axes[2].set_yticks([.001,.01,.1,1])
    axes[2].legend(loc='upper left',labelspacing=.4,handlelength=1.3,fontsize=8.5)
    for ax in axes[:2]: ax.set_xlim(0,100); ax.set_xticks([0,50,100])
    # Restore the active-window view in the approved original, not a new mask.
    axes[2].set_xlim(80,85); axes[2].set_xticks([80,82,84])

    # Full (F,G) projection; nearby stationary states are drawn in an explicitly
    # linked zoom below this panel, never displaced in data coordinates.
    branch=diag['branches']; bf=a(branch['F']); bg=a(branch['freeEnergyTotal']); bs=np.asarray(branch['stable'],bool)
    def projection(ax, arrows=False):
        for case,color,ls in zip(reps,colors,styles):
            ff,gg=a(case['F']),a(case['freeEnergy']['total'])
            ax.plot(ff,gg,color=COLORS[color],ls=ls,lw=1.3)
            ax.plot(ff[0],gg[0],'o',mfc='white',mec=COLORS[color],ms=3.5)
            ax.plot(ff[-1],gg[-1],'o',color=COLORS[color],ms=3.5)
            if arrows:
                length=np.hypot(np.diff(ff),np.diff(gg)); cumulative=np.cumsum(length)
                j=int(np.searchsorted(cumulative,.6*cumulative[-1])); k=min(j+2,len(ff)-1)
                ax.annotate('',xy=(ff[k],gg[k]),xytext=(ff[max(j-2,0)],gg[max(j-2,0)]),
                            arrowprops=dict(arrowstyle='->',color=COLORS[color],lw=.9))
        for i,(ff,gg) in enumerate(zip(bf,bg)):
            color='blue' if i==0 else 'red' if bs[i] else 'purple'
            ax.plot(ff,gg,marker='D' if bs[i] else 'x',color=COLORS[color],ms=5,ls='none',zorder=6)
    projection(axes[3], True)
    _finish(axes[3],r'$F$',r'$\mathcal{G}$')
    axes[3].set_xlim(-1.6,.43); axes[3].set_ylim(-.44,.058)
    axes[3].set_xticks([-1.5,-.5,.4])
    zoomd=fig.add_subplot(lower[1,0])
    projection(zoomd)
    _finish(zoomd,r'$F$',r'$\mathcal{G}$')
    zoomd.set_xlim(.02,.38); zoomd.set_ylim(.026,.036)
    zoomd.set_xticks([.05,.2,.35]); zoomd.set_yticks([.028,.034])
    zoomd.tick_params(labelsize=8.5)
    zoomd.text(.96,1.03,'detail',transform=zoomd.transAxes,ha='right',va='bottom',fontsize=8.5)
    indicator=axes[3].indicate_inset_zoom(zoomd,edgecolor='#777777',alpha=.5)
    for connector in indicator.connectors: connector.set_visible(False)
    open_start=_line('black','none','o'); open_start.set_markerfacecolor('white')
    filled_end=_line('black','none','o')
    axes[3].legend(handles=[(open_start,filled_end),
                        _line('blue','none','D','stable left'),_line('red','none','D','stable right'),
                        _line('purple','none','x','unstable')],
               labels=['start / end','stable left','stable right','unstable'],
               handler_map={tuple:HandlerTuple(ndivide=None)},
               loc='lower right',ncol=1,handlelength=1.3,fontsize=8,labelspacing=.35)

    p=np.array([c['P0'] for c in cases]); final=np.array([a(c['F'])[-1] for c in cases])
    coarse=a(diag['asymptoticBoundary']['intervals']); fine=a(diag['PbInterval'])
    def endpoint(ax):
        ax.plot(p,final,'o-',color=COLORS['black'],lw=1,ms=3.5)
        for i,ff in enumerate(bf):
            ax.axhline(ff,color=COLORS['blue' if i==0 else 'red' if bs[i] else 'purple'],
                       ls='--' if bs[i] else ':',lw=.9,zorder=0)
        ax.axvspan(*coarse,color=COLORS['gray'],alpha=.20,zorder=0)
        ax.axvline(coarse.mean(),color=COLORS['gray'],ls='-.',lw=.8)
        if not REMOVE_FITTED_FLUX_MARKER:
            ax.axvspan(*fine,color=COLORS['green'],alpha=.25,zorder=0)
            ax.axvline(fine.mean(),color=COLORS['green'],ls=':',lw=.9)
    endpoint(axes[4]); _finish(axes[4],r'$P_0$',r'$F(100)$')
    axes[4].set_xlim(-.002,.102); axes[4].set_xticks([0,.05,.1]); axes[4].set_ylim(-1.65,.32)
    zoome=fig.add_subplot(lower[1,1]); endpoint(zoome)
    _finish(zoome,r'$P_0$',r'$F(100)$')
    zoome.set_xlim(.07015,.07155); zoome.set_ylim(-1.65,.32)
    zoome.set_xticks([.0703,.0713]); zoome.set_xticklabels(['.0703','.0713'])
    zoome.set_yticks([-1.5,0]); indicator=axes[4].indicate_inset_zoom(zoome,edgecolor='#777777',alpha=.5)
    zoome.text(.96,1.03,'detail',transform=zoome.transAxes,ha='right',va='bottom',fontsize=8.5)
    for connector in indicator.connectors: connector.set_visible(False)
    legend_handles = [
        (Patch(fc=COLORS['gray'],alpha=.20), _line('gray','-.')),
        (Patch(fc=COLORS['green'],alpha=.25), _line('green',':'))]
    fig.legend(handles=legend_handles,
                   labels=['upwind basin bracket', 'fitted-flux basin bracket'],
                   handler_map={tuple:HandlerTuple(ndivide=1)},
                   loc='center',bbox_to_anchor=(.685,.649),ncol=2,
                   handlelength=1.8,columnspacing=1.4,fontsize=8)

    distance=np.array([c['minDistanceToUnstable'] for c in cases]); residence=np.array([c['unstableResidenceTime'] for c in cases])
    ftop.plot(p,distance,'o-',color=COLORS['blue'],ms=3.5)
    _finish(ftop,r'$P_0$','Minimum distance (0–100)')
    ftop.set_xlim(-.002,.102); ftop.set_xticks([0,.05,.1]); ftop.set_ylim(0,.67)
    ftop.yaxis.label.set_color(COLORS['blue']); ftop.tick_params(axis='y',colors=COLORS['blue'])
    ftop.text(.97,.97,'upwind',transform=ftop.transAxes,ha='right',va='top',fontsize=8.5)
    twin=ftop.twinx(); twin.plot(p,residence,'s--',color=COLORS['red'],ms=3.5,lw=1.1)
    clean_axis(twin); twin.spines['right'].set_visible(True); twin.spines['right'].set_color(COLORS['red'])
    twin.set_ylabel('Time within distance 0.1',color=COLORS['red'],labelpad=2)
    twin.tick_params(axis='y',colors=COLORS['red']); twin.set_ylim(0,80); twin.set_yticks([0,40,80])
    # Fine prediction uses its own recorded Pb and independent eigenvalue; the
    # distance column is consumed directly (no substitution of the upwind Pb).
    for side,color,ls in [('lower','blue','-'),('upper','orange','--')]:
        rows=sorted([r for r in predictions if r['side']==side],key=lambda r:float(r['distance']))
        xx=np.array([float(r['distance']) for r in rows]); yy=np.array([float(r['predictedResidence']) for r in rows])
        fbottom.plot(xx,yy,color=COLORS[color],ls=ls,lw=1.2)
        fbottom.fill_between(xx,[float(r['predictionLower']) for r in rows],
                            [float(r['predictionUpper']) for r in rows],color=COLORS[color],alpha=.14,lw=0)
        for role,marker in [('calibration','s'),('holdout','o')]:
            rr=[r for r in rows if r['role']==role]
            fbottom.plot([float(r['distance']) for r in rr],[float(r['measuredResidence']) for r in rr],
                         ls='none',marker=marker,color=COLORS[color],mfc=COLORS[color] if role=='calibration' else 'white',ms=4)
    _finish(fbottom,r'$\delta P_0$',r'$t_{\mathrm{res}}$')
    fbottom.set_xscale('log'); fbottom.set_xticks([1e-5,1e-4]); fbottom.set_ylim(30,92); fbottom.set_yticks([40,60,80])
    fbottom.xaxis.set_minor_locator(NullLocator())
    holdout_handle=_line('black','none','o','holdout'); holdout_handle.set_markerfacecolor('white')
    fbottom.legend(handles=[_line('blue','-',label='lower side'),_line('orange','--',label='upper side'),
                            _line('black','none','s','calibration'),holdout_handle],
                    loc='upper right',ncol=1,labelspacing=.25,handlelength=1.3,fontsize=8)
    for label,ax in zip('abcdef',axes): panel_label(ax,label)
    checks=dict(source_sha256=d['source_sha256'], prediction_sha256=hashlib.sha256(predfile.read_bytes()).hexdigest(),
                representative_P0=[c['P0'] for c in reps], all_P0_sorted=p.tolist(), coarse_boundary=coarse.tolist(),
                fine_boundary=fine.tolist(), dissipation_case_P0=critical['P0'], active_rate_samples=int(active.sum()),
                rate_display_time_window=[80,85],displayed_active_rate_samples=int((active&(tt>=80)&(tt<=85)).sum()),
                max_component_sum_residual=float(np.max(np.abs(a(free['total'])-a(free['entropy'])-a(free['localPotential'])-a(free['interaction'])))),
                max_representative_energy_increase=float(max(np.max(np.diff(a(c['freeEnergy']['total']))) for c in reps)),
                max_active_relative_rate_error=float(np.max(np.abs(minus[active]-rhs[active])/np.maximum(minus[active],rhs[active]))),
                median_active_relative_rate_error=float(np.median(np.abs(minus[active]-rhs[active])/np.maximum(minus[active],rhs[active]))),
                original_spatial_grid_points=len(diag['x']),
                original_rhs_rate_method=free['rhsRateMethod'],
                prediction_rows=len(predictions), prediction_band='Deterministic propagation of basin, eigenvalue and calibration-event intervals; not statistical confidence.',
                local_prediction_metadata=local_meta,
                projection_uses_absolute_G=True, final_curve_sorted_without_data_change=True)
    assert checks['max_component_sum_residual'] < 1e-12
    (folder/'fig5_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    (folder/'residence_prediction.csv').write_bytes(predfile.read_bytes())
    save_figure(fig,out,'fig5_free_energy_basin_selection')
    _manifest(root, folder, out)
    return checks


