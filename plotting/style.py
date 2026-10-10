"""Approved figure layout reused to render the fresh numerical results."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.legend import Legend

COLORS = dict(blue='#0072B2',orange='#D99500',gray='#696969',black='#171717',
              red='#D55E00',green='#009E73',purple='#8D579D',sky='#56B4E9')

def setup_style():
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
        'axes.titlesize':9,'legend.fontsize':9,'xtick.labelsize':8.5,'ytick.labelsize':8.5,
        'axes.linewidth':.8,'lines.linewidth':1.3,'lines.markersize':3.5,
        'xtick.direction':'out','ytick.direction':'out','xtick.major.width':.8,
        'ytick.major.width':.8,'xtick.major.size':3,'ytick.major.size':3,
        'axes.grid':False,'legend.frameon':True,'legend.framealpha':.30,
        'legend.facecolor':'white','legend.edgecolor':'#999999','legend.fancybox':False,
        'legend.borderpad':.35,'legend.labelspacing':.3,'legend.handletextpad':.5,
        'figure.facecolor':'white',
        'axes.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,
        'ps.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom',
        'mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold',
        'mathtext.cal':'STIXGeneral','mathtext.fallback':'stix',
        'figure.dpi':130,'savefig.dpi':300,'axes.unicode_minus':True})

def panel_label(ax,label):
    label=str(label).strip('()')
    # Position relative to the actual panel, not an exported bitmap.
    return ax.text(-.16,1.035,f'({label})',transform=ax.transAxes,
                   ha='left',va='bottom',fontsize=11,fontweight='bold',clip_on=False)

def clean_axis(ax,heatmap=False):
    ax.grid(False)
    ax.tick_params(direction='out',labelsize=8.5)
    ax.spines['top'].set_visible(heatmap)
    ax.spines['right'].set_visible(heatmap)
    for spine in ax.spines.values(): spine.set_linewidth(.8)

def save_figure(fig,out_dir,stem):
    out_dir=Path(out_dir).resolve()
    if any(part.lower()=='paper_text' for part in out_dir.parts):
        raise ValueError('Candidate exporter cannot write inside paper_text')
    out_dir.mkdir(parents=True,exist_ok=True)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    bounds=fig.bbox
    findings=[]
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text().strip(): continue
        box=item.get_window_extent(renderer)
        findings.append({'text':item.get_text(),'size_pt':item.get_fontsize(),
            'outside_canvas': bool(box.x0 < -1 or box.y0 < -1 or box.x1 > bounds.x1+1 or box.y1 > bounds.y1+1),
            'bbox_px':[float(v) for v in box.extents]})
    legend_records=[]
    for legend in fig.findobj(match=Legend):
        box=legend.get_window_extent(renderer)
        ax=legend.axes
        axbox=ax.get_window_extent(renderer) if ax is not None else None
        legend_records.append({'labels':[t.get_text() for t in legend.get_texts()],
            'title':legend.get_title().get_text(),'white_frame_alpha':legend.get_frame().get_alpha(),
            'bbox_px':[float(v) for v in box.extents],
            'axis_bbox_px':None if axbox is None else [float(v) for v in axbox.extents],
            'inside_axis':None if axbox is None else bool(box.x0>=axbox.x0-1 and box.y0>=axbox.y0-1 and box.x1<=axbox.x1+1 and box.y1<=axbox.y1+1)})
    paths=[]
    for ext in ('pdf','svg','png'):
        path=out_dir/f'{stem}.{ext}'
        fig.savefig(path,bbox_inches=None,dpi=300)
        paths.append({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    (out_dir/f'{stem}_layout.json').write_text(json.dumps({'width_mm':fig.get_figwidth()*25.4,
        'height_mm':fig.get_figheight()*25.4,'text':findings,'legends':legend_records,'outputs':paths},indent=2),encoding='utf-8')
    plt.close(fig)
