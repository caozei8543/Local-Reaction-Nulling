"""Frequency trajectories (Fig. 6, Fig. 9) from actual exported sweep rows."""
import itertools
from pathlib import Path
import numpy as np
from .core import *
from .io import *
from .plots import save
import matplotlib.pyplot as plt

def sweeps(path,out):
    out=Path(out);rs=rows(path,['design','source','eta_l','eta_r','frequency_hz','s11_re','s11_im','s21_re','s21_im']);groups={}
    for r in rs:
        k=(r['design'],r['source'],number(r['eta_l']),number(r['eta_r']));f=number(r['frequency_hz'])
        if k[0] not in DESIGNS or k[1] not in ('simulation','measurement') or f<=0: raise ValueError('Invalid sweep metadata')
        g=groups.setdefault(k,{})
        if f in g: raise ValueError(f'Duplicate frequency {k}/{f}')
        g[f]=(zrow(r,'s11'),zrow(r,'s21'))
    def curve(k):
        g=groups[k];f=np.array(sorted(g));v=np.array([g[x] for x in f]);return f,v
    tab=[];fig,axs=plt.subplots(1,3,figsize=(12,4),sharey=True)
    for ax,d in zip(axs,DESIGNS):
        keys=sorted(k for k in groups if k[0]==d and k[1]=='simulation' and k[3]==0 and -.1<=k[2]<=.1)
        if not keys or (d,'simulation',0.,0.) not in groups: raise ValueError(f'{d}: baseline and left-only sweep missing')
        f0,v0=curve((d,'simulation',0.,0.));mask0=(f0>=5.5e9)&(f0<=6.1e9)
        if not mask0.any(): raise ValueError('No samples in nominal sweep interval')
        ref=f0[mask0][np.argmin(abs(v0[mask0,1]))];shifts=[];x=[];z=[]
        for k in keys:
            f,v=curve(k);mask=(f>=5.5e9)&(f<=6.1e9);f,v=f[mask],v[mask];
            if len(f)<2 or not np.array_equal(f,f0[mask0]): raise ValueError('Fig. 6 requires common frequency grid')
            fm=f[np.argmin(abs(v[:,1]))];shift=(fm-ref)/1e6;shifts.append(abs(shift));x.append(k[2]);z.append(db(v[:,1]));tab.append(dict(design=d,eta_l=k[2],fmin_hz=fm,shift_mhz=shift))
        if len(x)<2: raise ValueError('Need more than one load')
        im=ax.pcolormesh(x,f/1e9,np.array(z).T,shading='auto',vmin=-30,vmax=-10);ax.set_title(f'{d}\nmax sampled shift {max(shifts):.2f} MHz');ax.set_xlabel('Left perturbation');ax.set_ylim(5.5,6.1)
    axs[0].set_ylabel('Frequency (GHz)');fig.colorbar(im,ax=list(axs),label='S21 (dB)');save(fig,out/'Fig6_from_sweeps.png');write_csv(out/'fmin_shifts.csv',tab)
    k=('RN','measurement',0.,0.)
    if k in groups:
        fig,ax=plt.subplots(figsize=(7,4))
        for typ,ls in [('simulation','-'),('measurement','--')]:
            f,v=curve(('RN',typ,0.,0.));ax.plot(f/1e9,db(v[:,0]),ls,label=f'{typ} S11');ax.plot(f/1e9,db(v[:,1]),ls,label=f'{typ} S21')
        ax.set(xlabel='Frequency (GHz)',ylabel='Magnitude (dB)',title='RN nominal supplied sweeps');ax.grid(alpha=.2);ax.legend();save(fig,out/'Fig9_from_sweeps.png')
    return tab
