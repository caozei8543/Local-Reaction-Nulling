"""Publication-oriented reconstructions, explicitly labelled by provenance."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .core import DESIGNS,CONDITIONS,db

def save(fig,path):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=220,bbox_inches='tight');plt.close(fig)

def reaction_plot(a,b,path,label='Reconstruction from reported integrals'):
    fig,ax=plt.subplots(figsize=(6,5));ax.axhline(0,color='.65',lw=.8);ax.axvline(0,color='.65',lw=.8)
    for z,color,name in [(a,'tab:red','Left subregion a'),(b,'tab:blue','Left subregion b'),(a+b,'tab:green','Residual')]:
        ax.annotate('',(z.real,z.imag),(0,0),arrowprops=dict(arrowstyle='->',color=color,lw=2));ax.plot([],[],color=color,label=f'{name}: {abs(z):.4f}')
    ax.plot([a.real,(a+b).real],[a.imag,(a+b).imag],'--',color='tab:blue',alpha=.7)
    ax.scatter([(a+b).real],[(a+b).imag],color='tab:green',zorder=3)
    ax.set(xlabel='Real part (ND normalized)',ylabel='Imaginary part (ND normalized)',title=label);ax.set_aspect('equal',adjustable='datalim');ax.legend(loc='best',fontsize=8);ax.grid(alpha=.2);save(fig,path)

def comparison_plot(d,path,errors=None,label='Reported ratios (no raw error bars)'):
    fig,ax=plt.subplots(figsize=(8,4.5));metrics=['FD left','FD right','Common load','Unequal load'];x=np.arange(4)
    for i,k in enumerate(DESIGNS):
        v=d[k];y=[v['fd_left_mag_deg'][0],v['fd_right_mag_deg'][0],*[v['meas_finite_ratios'][q] for q in CONDITIONS]]
        ax.scatter(x+(i-1)*.16,y,label=k,marker=['s','^','o'][i],s=54)
        if errors is not None: ax.errorbar(x[2:]+(i-1)*.16,y[2:],yerr=errors[k],fmt='none',capsize=3)
    ax.set_xticks(x,metrics);ax.set_yscale('log');ax.set_ylabel('Normalized magnitude');ax.set_title(label);ax.grid(axis='y',which='both',alpha=.25);ax.legend();save(fig,path)

def height_integrated_density_plot(density,path):
    """Fig.3-style height integration for a uniform midpoint tensor grid; no spatial interpolation."""
    data=[r for r in density if r['design']=='RN'];xy={};areas={}
    for region in ('La','Lb','Ra','Rb'):
        rs=[r for r in data if r['region']==region];xs=sorted({r['x_m'] for r in rs});ys=sorted({r['y_m'] for r in rs});zs=sorted({r['z_m'] for r in rs})
        if len(rs)!=len(xs)*len(ys)*len(zs):return False
        # Known paper height = 2 mm, tensor-product uniform grid required.
        w=np.array([r['weight_m3'] for r in rs]);
        if not np.allclose(w,w[0],rtol=1e-6,atol=0):return False
        dz=.002/len(zs)
        for r in rs:
            key=(r['x_m'],r['y_m']);xy[key]=xy.get(key,0)+r['reaction_density_re']*dz
    if not xy:return False
    a=np.array([(x*1e3,y*1e3,v) for (x,y),v in xy.items()]);vmax=max(abs(a[:,2]))
    if vmax==0:return False
    fig,ax=plt.subplots(figsize=(7,4));im=ax.scatter(a[:,0],a[:,1],c=a[:,2]/vmax,cmap='RdBu_r',vmin=-1,vmax=1,marker='s',s=25);fig.colorbar(im,ax=ax,label='Normalized height-integrated Re reaction density');ax.set(xlabel='x (mm)',ylabel='y (mm)',title='RN supplied tensor-grid fields');ax.set_aspect('equal');save(fig,path);return True
