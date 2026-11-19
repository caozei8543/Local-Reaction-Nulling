"""Read single-ended 2-port Touchstone v1, 50-ohm S, RI/MA/DB (S11,S21,S12,S22)."""
from pathlib import Path
import numpy as np
from .core import F0,polar
from .io import rows,write_csv,zcols

def read_s2p(path):
    opts=None;tokens=[]
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line=line.split('!',1)[0].strip()
        if not line:continue
        if line.startswith('['):raise ValueError('Touchstone v2 not supported; export v1 .s2p')
        if line.startswith('#'):
            if opts is not None or tokens:raise ValueError('Unexpected option line')
            opts=line[1:].upper().split();continue
        tokens.extend(line.split())
    if opts is None:raise ValueError('Explicit Touchstone option line required')
    if len(opts)!=5 or opts[1]!='S' or opts[2] not in ('RI','MA','DB') or opts[3]!='R' or float(opts[4])!=50:raise ValueError('Require 50-ohm single-ended S data')
    scale={'HZ':1,'KHZ':1e3,'MHZ':1e6,'GHZ':1e9}[opts[0]]
    if not tokens or len(tokens)%9:raise ValueError('Invalid two-port Touchstone record length')
    a=np.array(tokens,float).reshape(-1,9);f=a[:,0]*scale
    if not np.isfinite(a).all() or np.any(np.diff(f)<=0):raise ValueError('Nonfinite or unsorted frequency records')
    x,y=a[:,1::2],a[:,2::2]
    s=x+1j*y if opts[2]=='RI' else (x if opts[2]=='MA' else 10**(x/20))*np.exp(1j*np.deg2rad(y))
    return f,s

def s21_at_f0(path):
    f,s=read_s2p(path);i=np.where(abs(f-F0)<=1)[0]
    if len(i)!=1:raise ValueError(f'{path}: exactly one 5.8 GHz sample required (no interpolation)')
    return complex(s[i[0],1])

def convert_manifest(manifest,output):
    manifest=Path(manifest);rs=rows(manifest,['design','condition','cycle','loaded_s2p','unloaded_s2p']);out=[]
    for r in rs:
        a=s21_at_f0(manifest.parent/r['loaded_s2p']);b=s21_at_f0(manifest.parent/r['unloaded_s2p']);out.append(dict(design=r['design'],condition=r['condition'],cycle=r['cycle'],frequency_hz=F0,**zcols('loaded',a),**zcols('unloaded',b)))
    write_csv(output,out);return out
