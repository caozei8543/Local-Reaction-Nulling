"""Generate common quadrature coordinates for two HFSS port-field exports.
Geometry locations MUST be supplied from the original CAD, not guessed from the figure.
"""
import itertools
import numpy as np
from .io import read_json,write_csv
from .raw import REGIONS

def generate(config,output):
    m=read_json(config);regions=m['regions_mm'];result=[]
    for name in REGIONS:
        r=regions[name];lo=np.array(r['min'],float);hi=np.array(r['max'],float);n=np.array(r['cells'],int)
        if lo.shape!=(3,) or hi.shape!=(3,) or n.shape!=(3,) or np.any(hi<=lo) or np.any(n<1) or not np.allclose(sorted(hi-lo),[2,3,9],rtol=0,atol=1e-8): raise ValueError(f'{name}: supply actual 9×3×2 mm CAD box and positive cell counts')
        step=(hi-lo)/n;axes=[lo[i]+(np.arange(n[i])+.5)*step[i] for i in range(3)];w=float(np.prod(step))*1e-9
        for x,y,z in itertools.product(*axes):result.append(dict(region=name,x_m=x*1e-3,y_m=y*1e-3,z_m=z*1e-3,weight_m3=w,g=1))
    write_csv(output,result);return result

def merge(grid_path,port1,port2,design,output):
    """Merge exported Cartesian fields by exact point_id from user-aligned exports."""
    from .io import rows,unique_map,number,zrow,zcols
    from .core import DESIGNS,F0
    if design not in DESIGNS:raise ValueError('Unknown design')
    grid=rows(grid_path,['region','x_m','y_m','z_m','weight_m3','g']);key=lambda r:tuple(number(r[k]) for k in ('x_m','y_m','z_m'))
    required=['x_m','y_m','z_m']+[f'E{a}_{c}' for a in 'xyz' for c in ('re','im')];a,b=[unique_map(rows(p,required),key) for p in (port1,port2)]
    ks={key(r) for r in grid}
    if len(ks)!=len(grid) or set(a)!=ks or set(b)!=ks:raise ValueError('Both port exports must have exactly the same coordinates as the grid, in metres')
    result=[]
    for r in grid:
        k=key(r);line=dict(design=design,frequency_hz=F0,**r)
        for p,rec in enumerate((a[k],b[k]),1):
            for ax in 'xyz':line.update(zcols(f'E{p}{ax}',zrow(rec,f'E{ax}')))
        result.append(line)
    write_csv(output,result);return result
