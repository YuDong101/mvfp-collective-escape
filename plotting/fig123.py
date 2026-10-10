"""Frozen and dynamic candidate figures from the dedicated MATLAB export."""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from .style import setup_style, panel_label, clean_axis, save_figure, COLORS

STEMS={1:'fig1_frozen_landscape_limitations',2:'fig2_dynamic_feedback_mechanism',3:'fig3_finite_time_readout_transition'}
P0S=[0,.04,.06,.0647,.07,.08,.10,.16]
P0COLORS=['#440154','#46327E','#365C8D','#277F8E','#1F9E89','#4CB74C','#80AF27','#B99400']
STYLES=['-','--','-.',':','-','--','-.',':']

def color(p):
    return P0COLORS[int(np.argmin(np.abs(np.array(P0S)-float(p))))]

def rows(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        result=list(csv.DictReader(f))
    for r in result:
        for k,v in r.items():
            try: r[k]=float(v)
            except (ValueError,TypeError): pass
    return result

def arr(rs,key): return np.asarray([r[key] for r in rs],dtype=float)
def select(rs,key,val): return [r for r in rs if abs(float(r[key])-val)<1e-10]
def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def number_array(value): return np.asarray(value,dtype=float)

def data_dir(root): return Path(__file__).resolve().parents[1]/'data/fig123'

def provenance(folder,n,files,notes):
    return {'figure':n,'notes':notes,'sources':[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]}


def top_colorbar(fig,ax,im,title,ticks=None,extend='neither'):
    b=ax.get_position()
    cax=fig.add_axes([b.x0,b.y1+.035,b.width,.012])
    cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=ticks,extend=extend)
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.tick_params(labelsize=8,pad=1,length=2)
    cb.ax.set_title(title,fontsize=9,pad=13)
    return cb

def build_figure1(root:Path,out:Path):
    setup_style();folder=data_dir(root)
    summary=rows(folder/'fig1_frozen_summary.csv');summary.sort(key=lambda r:r['P0'])
    pot=read(folder/'fig1_potentials.json')
    threshold=rows(folder/'fig3_threshold_brackets.csv')
    thr=[r for r in threshold if r['T']==100]
    interval=[min(r['P0Lower'] for r in thr),max(r['P0Upper'] for r in thr)]
    dynamic=sum(interval)/2
    fig=plt.figure(figsize=(179/25.4,133/25.4))
    axes=[fig.add_axes(b) for b in [[.11,.57,.31,.31],[.62,.57,.31,.31],[.11,.13,.31,.31],[.62,.13,.31,.31]]]
    x=number_array(pot['x']);us=number_array(pot['UminusMin']);ps=np.ravel(number_array(pot['P0']))
    palette=['#440154','#365C8D','#21918C','#61B64A','#B99400']
    handles=[]
    for k,p in enumerate(ps):
        h,=axes[0].plot(x,us[:,k],color=palette[k],ls=['-','--','-.',':',(0,(5,1,1,1,1,1))][k],label=fr'$P_0={p:g}$')
        handles.append(h)
    axes[0].set(xlim=(-1,1),ylim=(0,1),xlabel='$x$',ylabel=r'$V_0$')
    axes[0].legend(handles=handles,ncol=2,loc='upper center',fontsize=8,columnspacing=.7,handlelength=1.5)
    p=arr(summary,'P0')
    ext=[('xLeftMinimum','left min.',COLORS['blue'],'-'),('xSaddle','saddle',COLORS['gray'],'--'),('xRightMinimum','right min.',COLORS['orange'],'-'),('xSingleMinimum','single min.',COLORS['green'],'-.')]
    for key,label,c,ls in ext: axes[1].plot(p,arr(summary,key),color=c,ls=ls,label=label)
    axes[1].set(xlim=(0,1),xlabel='$P_0$',ylabel='extremum $x$')
    axes[1].legend(loc='center right',fontsize=8,ncol=1,handlelength=1.7)
    twins=[]
    for i,lk,rk,ylabel in [(2,'barrierRightToLeft','barrierLeftToRight',r'$B_{R\to L}$'),(3,'log10TauRightToSaddle','log10TauLeftToSaddle',r'$L_{R\to L}$')]:
        ax=axes[i];twin=ax.twinx();twins.append(twin)
        left=arr(summary,lk);right=arr(summary,rk)
        if i==2: left=left/.005;right=right/.005
        else: left=left-2;right=right-2
        ax.plot(p,left,color=COLORS['blue']);twin.plot(p,right,color=COLORS['orange'],ls='--')
        maxp=float(p[np.isfinite(left)&np.isfinite(right)].max())
        ax.set(xlim=(0,maxp),xlabel='$P_0$',ylabel=ylabel)
        twin.set_ylabel(r'$B_{L\to R}$' if i==2 else r'$L_{L\to R}$',color=COLORS['orange'])
        ax.yaxis.label.set_color(COLORS['blue']);ax.tick_params(axis='y',colors=COLORS['blue'])
        twin.tick_params(axis='y',colors=COLORS['orange']);twin.spines['top'].set_visible(False)
        twin.spines['left'].set_visible(False);twin.spines['right'].set_color(COLORS['orange'])
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4));twin.yaxis.set_major_locator(MaxNLocator(4))
    tau=arr(summary,'log10TauRightToSaddle')-2
    finite=np.isfinite(tau);px=p[finite];tv=tau[finite]
    cross=np.where(tv[:-1]*tv[1:]<=0)[0]
    frozen=None
    if len(cross):
        frozen=float(read(folder/'fig1_metadata.json')['frozenTauRLEquals100P0'])
    axes[3].axhline(0,color=COLORS['gray'],ls=':',lw=.9)
    axes[3].legend(handles=[Line2D([],[],color=COLORS['gray'],ls=':',label='left-axis zero')],loc='upper center',fontsize=8,handlelength=1.5)
    for label,ax in zip('abcd',axes):clean_axis(ax);panel_label(ax,label)
    save_figure(fig,out,STEMS[1])
    return provenance(folder,1,[folder/f for f in ['fig1_frozen_summary.csv','fig1_potentials.json','fig1_metadata.json','fig3_threshold_brackets.csv']],{'dynamic_interval':interval,'frozen_crossing_refined_root':frozen,'D':.005,'T':100,'no_band_inflation':True})

def build_figure2(root:Path,out:Path):
    setup_style();folder=data_dir(root)
    full=rows(folder/'fig2_dynamic_full_traces.csv');controls=rows(folder/'fig2_same_density_controls.csv')
    events=rows(folder/'fig2_event_timing.csv');pot=read(folder/'fig2_dynamic_potential.json')
    fig=plt.figure(figsize=(179/25.4,125/25.4))
    axes=[fig.add_axes([left,bottom,.205,.30]) for bottom in [.62,.11] for left in [.105,.43,.755]]
    handles=[]
    for k,p in enumerate(P0S):
        rs=select(full,'P0',p);rs.sort(key=lambda r:r['t']);t=arr(rs,'t')
        for j,key in enumerate(['PFixed','F','barrierRightToLeft']):
            y=arr(rs,key)/(.005 if j==2 else 1)
            h,=axes[j].plot(t,y,color=P0COLORS[k],ls=STYLES[k],lw=1.3)
            if j==0:handles.append(Line2D([],[],color=P0COLORS[k],ls=STYLES[k],label=f'{p:g}'))
    axes[0].legend(handles=handles,title='$P_0$',ncol=2,loc='upper left',fontsize=8,title_fontsize=9,columnspacing=.65,handlelength=1.4)
    for j in range(3):axes[j].set(xlim=(0,100),xlabel='$t$');axes[j].set_xticks([0,50,100])
    axes[0].set(ylim=(-.02,1.04),ylabel=r'$P_{\mathrm{fixed}}(t)$')
    axes[1].set_ylabel('$F(t)$');axes[2].set_ylabel(r'$B_{R\to L}$')
    axes[0].set_yticks([0,.5,1]);axes[1].yaxis.set_major_locator(MaxNLocator(4));axes[2].yaxis.set_major_locator(MaxNLocator(4))
    t=number_array(pot['t']);x=number_array(pot['x']);u=number_array(pot['UminusMin'])
    im=axes[3].pcolormesh(t,x,u,cmap='magma',vmin=0,vmax=1,shading='nearest',rasterized=True)
    for key,col,ls in [('xLeftMinimum',COLORS['sky'],'-'),('xSaddle',COLORS['black'],'--'),('xRightMinimum',COLORS['green'],'-.'),('xSingleMinimum',COLORS['sky'],':')]:
        axes[3].plot(t,number_array(pot[key]),color=col,ls=ls,lw=1.1)
    axes[3].set(xlim=(0,100),ylim=(-1,1),xlabel='$t$',ylabel='$x$');axes[3].set_xticks([0,50,100]);axes[3].set_yticks([-1,0,1])
    assert float(pot['P0'])==.08
    top_colorbar(fig,axes[3],im,r'$V_t$',ticks=[0,.5,1],extend='max')
    axes[3].text(.96,.96,'$P_0=0.08$',transform=axes[3].transAxes,fontsize=8,va='top',ha='right',color='black')
    case_event=[r for r in events if abs(r['P0']-.08)<1e-10 and r['q']==.5][0]
    for p in [.08,.10,.16]:
        ev=select(events,'P0',p)[0]
        axes[2].plot((ev['lastBistableLower']+ev['lastBistableUpper'])/2,0,'x',color=color(p),ms=5,clip_on=False)
    levcolors=[COLORS['purple'],COLORS['blue'],COLORS['green'],COLORS['orange']]
    censored=set();eh=[]
    for q,c,m in zip([.2,.4,.6,.8],levcolors,['o','s','D','v']):
        rs=select(events,'q',q);rs.sort(key=lambda r:r['P0'])
        xq=arr(rs,'P0');y=arr(rs,'time');observed=np.isfinite(y)&(y<=100)
        axes[4].plot(xq,np.where(observed,y,np.nan),color=c,marker=m,lw=1,ms=3.5)
        censored.update(xq[~observed].tolist());eh.append(Line2D([],[],color=c,marker=m,label=f'$P={q:g}$',lw=1))
    foldrows=[]
    for p in P0S:
        rs=select(events,'P0',p)
        if not rs:continue
        lo=float(rs[0]['lastBistableLower']);hi=float(rs[0]['lastBistableUpper'])
        if np.isfinite(lo) and np.isfinite(hi) and hi<=100:
            mid=(lo+hi)/2;axes[4].errorbar(p,mid,yerr=(hi-lo)/2,color=COLORS['black'],fmt='s',ms=3,capsize=2)
            foldrows.append((p,mid))
        else:censored.add(p)
    if foldrows:axes[4].plot(*np.array(foldrows).T,color=COLORS['black'],ls='--',lw=.9)
    axes[4].plot(sorted(censored),np.full(len(censored),100),'^',mfc='white',mec=COLORS['black'],ms=4,ls='none',clip_on=False)
    eh += [Line2D([],[],color=COLORS['black'],ls='--',marker='s',label='well loss'),Line2D([],[],color=COLORS['black'],ls='none',marker='^',mfc='white',label='censored')]
    axes[4].legend(handles=eh,ncol=1,loc='lower left',fontsize=8,handlelength=1.2)
    axes[4].set(xlim=(-.004,.166),ylim=(0,106),xlabel='$P_0$',ylabel='event time');axes[4].set_xticks([0,.08,.16]);axes[4].set_yticks([0,50,100])
    ch=[]
    for key,c,ls,label in [('full',COLORS['blue'],'-','full'),('frozenFeedback',COLORS['orange'],'--','frozen'),('noCoupling',COLORS['gray'],':','uncoupled')]:
        rs=[r for r in controls if abs(r['P0']-.08)<1e-10 and r['control']==key];rs.sort(key=lambda r:r['t'])
        h,=axes[5].plot(arr(rs,'t'),arr(rs,'PFixed'),color=c,ls=ls,label=label);ch.append(h)
    axes[5].set(xlim=(0,100),ylim=(-.02,1.05),xlabel='$t$',ylabel=r'$P_{\mathrm{fixed}}(t)$');axes[5].set_xticks([0,50,100]);axes[5].set_yticks([0,.5,1])
    axes[5].axhline(.5,color=COLORS['gray'],ls=':',lw=.8)
    axes[5].text(.98,.04,'$P_0=0.08$',transform=axes[5].transAxes,ha='right',fontsize=8)
    inset=axes[5].inset_axes([.20,.63,.50,.25])
    for j,(key,c) in enumerate([('full',COLORS['blue']),('frozenFeedback',COLORS['orange']),('noCoupling',COLORS['gray'])]):
        rs=[r for r in controls if abs(r['P0']-.06)<1e-10 and r['control']==key];rs.sort(key=lambda r:r['t'])
        inset.plot(j,rs[-1]['PFixed'],'o',color=c,ms=3.5)
    inset.set(xlim=(-.5,2.5),ylim=(0,.25),xticks=[],yticks=[0,.2]);inset.tick_params(labelsize=8,length=2,pad=1)
    inset.set_title('$P_0=0.06$; $P(T)$',fontsize=8,pad=1);clean_axis(inset)
    axes[5].legend(handles=ch,loc='center left',bbox_to_anchor=(.02,.37),fontsize=8,handlelength=1.5)
    ext=[Line2D([],[],color=c,ls=ls,label=lab) for c,ls,lab in [(COLORS['sky'],'-','left min.'),(COLORS['black'],'--','saddle'),(COLORS['green'],'-.','right min.'),(COLORS['sky'],':','single min.')]]
    axes[3].legend(handles=ext,ncol=1,loc='upper left',fontsize=8,handlelength=1.5,labelspacing=.25)
    for i,ax in enumerate(axes):clean_axis(ax,heatmap=i==3);panel_label(ax,'abcdef'[i])
    save_figure(fig,out,STEMS[2])
    return provenance(folder,2,[folder/f for f in ['fig2_dynamic_full_traces.csv','fig2_same_density_controls.csv','fig2_event_timing.csv','fig2_dynamic_potential.json']],{'P0':P0S,'potential_representative':.08,'control_representative':.08,'counterexample':.06,'censoring':'triangle at horizon denotes one or more unobserved displayed events; absent colored curves identify missing levels','potential_saturation':1})

def build_figure3(root:Path,out:Path):
    setup_style();folder=data_dir(root)
    initial=rows(folder/'fig3_initial_densities.csv');rep=rows(folder/'fig3_representative_readout_traces.csv')
    thr=rows(folder/'fig3_threshold_brackets.csv');cens=rows(folder/'fig3_t50_censoring_T100.csv')
    density=read(folder/'fig3_density_evolution.json');heat=read(folder/'fig3_readout_heatmap.json')
    fig=plt.figure(figsize=(179/25.4,135/25.4))
    axes=[fig.add_axes([left,bottom,.205,.27]) for bottom in [.565,.16] for left in [.105,.43,.755]]
    hs=[]
    for p,ls in zip([.06,.07,.08],['-','--','-.']):
        rs=select(initial,'P0',p);rs.sort(key=lambda r:r['x'])
        h,=axes[0].plot(arr(rs,'x'),arr(rs,'rho0'),color=color(p),ls=ls,label=fr'$P_0={p:g}$');hs.append(h)
        rs=select(rep,'P0',p);rs.sort(key=lambda r:r['t']);t=arr(rs,'t');y=arr(rs,'PFixed')
        axes[3].plot(t,y,color=color(p),ls=ls)
    axes[0].set(xlim=(-.65,.65),xlabel='$x$',ylabel=r'$\rho_0(x)$');axes[0].set_xticks([-.5,0,.5]);axes[0].yaxis.set_major_locator(MaxNLocator(3))
    axes[3].legend(handles=hs,labels=['0.06','0.07','0.08'],title='$P_0$',ncol=1,loc='upper left',fontsize=8.5,handlelength=1.5)
    zax=axes[0].inset_axes([.28,.48,.42,.40])
    for p,ls in zip([.06,.07,.08],['-','--','-.']):
        rs=select(initial,'P0',p);rs.sort(key=lambda r:r['x']);zax.plot(arr(rs,'x'),arr(rs,'rho0'),color=color(p),ls=ls,lw=1)
    zax.set(xlim=(-.54,-.40),ylim=(0,1.8),xticks=[-.5],yticks=[0,1.5]);clean_axis(zax)
    axes[0].indicate_inset_zoom(zax,edgecolor=COLORS['gray'],alpha=.6)
    x=number_array(density['x']);t=number_array(density['t']);rho=number_array(density['rhoRaw'])
    if np.nanmin(rho)<-1e-6: raise ValueError('Raw negative density exceeds display offset')
    im=axes[1].pcolormesh(t,x,np.log10(rho+1e-6),cmap='viridis',vmin=-4,vmax=1.4,shading='nearest',rasterized=True)
    axes[1].set(xlim=(0,200),ylim=(-.65,.65),xlabel='$t$',ylabel='$x$');axes[1].set_xticks([0,100,200]);axes[1].set_yticks([-.5,0,.5])
    axes[1].axvline(100,color='white',ls=':',lw=.9)
    axes[1].text(.98,.97,'$P_0=0.07$',transform=axes[1].transAxes,ha='right',va='top',fontsize=8,color='white')
    top_colorbar(fig,axes[1],im,r'$L_\rho$',ticks=[-4,-2,0,1.4],extend='both')
    hp=number_array(heat['P0']);ht=number_array(heat['t']);values=number_array(heat['PFixed'])
    assert values.shape==(hp.size,ht.size)
    im=axes[2].pcolormesh(ht,hp,values,cmap='viridis',vmin=0,vmax=1,shading='nearest',rasterized=True)
    axes[2].contour(ht,hp,values,levels=[.5],colors='white',linewidths=.9)
    overlay=[]
    for T in [50,100]:
        rs=[r for r in thr if r['T']==T];lo=min(r['P0Lower'] for r in rs);hi=max(r['P0Upper'] for r in rs);mid=(lo+hi)/2
        axes[2].errorbar(T,mid,yerr=[[mid-lo],[hi-mid]],fmt='o',color=COLORS['red'],mfc='none',ms=4,capsize=2,clip_on=False)
        overlay.append({'T':T,'P0_interval':[lo,hi]})
    axes[2].set(xlim=(0,100),ylim=(0,.3),xlabel='$t$',ylabel='$P_0$');axes[2].set_xticks([0,50,100]);axes[2].set_yticks([0,.1,.2,.3])
    top_colorbar(fig,axes[2],im,r'$P_{\mathrm{fixed}}$',ticks=[0,.5,1])
    axes[3].set(xlim=(0,200),ylim=(-.02,1.04),xlabel='$t$',ylabel=r'$P_{\mathrm{fixed}}(t)$');axes[3].set_xticks([0,100,200]);axes[3].set_yticks([0,.5,1])
    axes[3].axhline(.5,color=COLORS['gray'],ls=':',lw=.9);axes[3].axvline(100,color=COLORS['gray'],ls=':',lw=.9)
    th=[]
    for nx,c,m in [(2001,COLORS['blue'],'o'),(4001,COLORS['orange'],'s')]:
        rs=select(thr,'Nx',nx);rs=[r for r in rs if np.isfinite(r['P0Lower'])];rs.sort(key=lambda r:r['T'])
        tt=arr(rs,'T');lo=arr(rs,'P0Lower');hi=arr(rs,'P0Upper');mid=(lo+hi)/2
        axes[4].errorbar(tt,mid,yerr=[mid-lo,hi-mid],fmt=m+'-',color=c,ms=4,mfc='white' if nx==4001 else c,lw=1,capsize=2)
        th.append(Line2D([],[],color=c,marker=m,mfc='white' if nx==4001 else c,label=f'$N_x={nx}$'))
    axes[4].set(xlim=(35,410),ylim=(0,.13),xlabel='$T$',ylabel=r'$P_{0c}$');axes[4].set_xticks([50,200,400]);axes[4].set_yticks([0,.06,.12])
    zax=axes[4].inset_axes([.50,.60,.44,.31])
    for j,(nx,c,m) in enumerate([(2001,COLORS['blue'],'o'),(4001,COLORS['orange'],'s')]):
        rs=[r for r in thr if r['Nx']==nx and r['T']==100];r=rs[0];lo=r['P0Lower'];hi=r['P0Upper'];mid=(lo+hi)/2
        zax.errorbar(j,mid,yerr=[[mid-lo],[hi-mid]],fmt=m,color=c,ms=3,capsize=2)
    zax.set(xlim=(-.6,1.6),xticks=[0,1],xticklabels=['2k','4k']);zax.set_yticks([.07045,.07055]);zax.tick_params(labelsize=8,pad=1,length=2)
    zax.set_title('$T=100$',fontsize=8,pad=2);zax.ticklabel_format(axis='y',style='plain',useOffset=False);clean_axis(zax)
    # Full-precision intervals remain data; the zoom alone magnifies the display scale.
    axes[4].text(.99,.18,'$T=400$:\nno crossing',transform=axes[4].transAxes,ha='right',va='bottom',fontsize=8)
    axes[4].legend(handles=th,ncol=1,loc='lower left',fontsize=8,handlelength=1.3)
    cens.sort(key=lambda r:r['P0']);p=arr(cens,'P0');t50=arr(cens,'t50Actual')
    observed=np.isfinite(t50)&(t50<=100)
    axes[5].plot(p,np.where(observed,t50,np.nan),color=COLORS['blue'],lw=1)
    axes[5].plot(p[observed],t50[observed],'o',color=COLORS['blue'],ms=3)
    axes[5].plot(p[~observed],np.full((~observed).sum(),100),'^',mfc='white',mec=COLORS['orange'],ms=3.5,clip_on=False)
    axes[5].axhline(100,color=COLORS['gray'],lw=.8,ls=':')
    axes[5].set(xlim=(-.005,.305),ylim=(0,106),xlabel='$P_0$',ylabel='$t_{50}$');axes[5].set_xticks([0,.15,.3]);axes[5].set_yticks([0,50,100])
    axes[5].legend(handles=[Line2D([],[],color=COLORS['blue'],marker='o',lw=1,label='reached'),Line2D([],[],color=COLORS['orange'],marker='^',mfc='white',ls='none',label='censored')],loc='upper right',bbox_to_anchor=(1,.78),fontsize=8,handlelength=1)
    for i,ax in enumerate(axes):clean_axis(ax,heatmap=i in [1,2]);panel_label(ax,'abcdef'[i])
    save_figure(fig,out,STEMS[3])
    return provenance(folder,3,[folder/f for f in ['fig3_initial_densities.csv','fig3_representative_readout_traces.csv','fig3_threshold_brackets.csv','fig3_t50_censoring_T100.csv','fig3_density_evolution.json','fig3_readout_heatmap.json']],{'overlay_data_coordinates':overlay,'heatmap_range':[0,100,0,.3],'density_time_range':[0,200],'density_transform':'log10(rhoRaw+1e-6), display [-4,1.4] with marked saturation','events_from_saved_numeric_traces':True})
