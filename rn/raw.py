"""Raw exported-data analysis. No fallback to reconstructed numbers."""
from pathlib import Path
import itertools
import numpy as np
from .core import *
from .io import *
from .plots import save

FIELD_KEYS=['design','frequency_hz','region','x_m','y_m','z_m','weight_m3','g']+[f'E{p}{a}_{c}' for p in (1,2) for a in 'xyz' for c in ('re','im')]
REGIONS=('La','Lb','Ra','Rb')
FIELD_VOLUME=9*3*2*1e-9

def fields(path,metadata,out):
    """Each row is one matched quadrature sample; both port fields must be co-located."""
    m=read_json(metadata); out=Path(out)
    if m['time_convention']!='exp(+j*w*t)' or m['reference_ohm']!=50 or m['excitation_type']!='incident_power' or m['unit_phase'] is not True or m['matched_other_port'] is not True:
        raise ValueError('Field metadata must certify article port/time conventions')
    if m['coordinate_units']!='m' or m['field_units']!='V/m' or m['weight_units']!='m^3': raise ValueError('SI units required')
    rs=rows(path,FIELD_KEYS); unique_map(rs,lambda r:(r['design'],r['region'],*(number(r[x]) for x in ('x_m','y_m','z_m'))))
    if any(r['design'] not in DESIGNS or r['region'] not in REGIONS for r in rs): raise ValueError('Unknown design/region')
    result={};tab=[];density=[]
    for k in DESIGNS:
        reg={}
        for q in REGIONS:
            data=[r for r in rs if r['design']==k and r['region']==q]
            if not data: raise ValueError(f'Missing {k}/{q} field rows')
            for r in data: check_frequency(r)
            w=np.array([number(r['weight_m3']) for r in data]); g=np.array([number(r['g']) for r in data]);
            if not np.allclose(g,1,atol=1e-12,rtol=0): raise ValueError('The prescribed uniform window profiles require g=1')
            if not np.isclose(w.sum(),FIELD_VOLUME,rtol=.01,atol=0): raise ValueError(f'{k}/{q}: volume weights must sum to 54 mm^3 within 1%')
            e=[np.array([[zrow(r,f'E{p}{a}') for a in 'xyz'] for r in data]) for p in (1,2)]
            rr=reaction(*e,w,g);reg[q]=rr
            den=np.sum(e[0]*e[1],axis=1)*g
            for r,dv in zip(data,den): density.append(dict(design=k,region=q,x_m=number(r['x_m']),y_m=number(r['y_m']),z_m=number(r['z_m']),weight_m3=number(r['weight_m3']),**zcols('reaction_density',dv)))
        reg.update(L=reg['La']+reg['Lb'],R=reg['Ra']+reg['Rb']);result[k]=reg
    ref=result['ND']['L'];
    if ref==0: raise ValueError('ND reference reaction is zero')
    for k,reg in result.items():
        for q,r in reg.items():
            d=derivative_from_reaction(r,phasor=m['phasor'],p1_w=number(m['port1_incident_power_w']),p2_w=number(m['port2_incident_power_w']))
            tab.append(dict(design=k,region=q,**zcols('reaction',r),**zcols('normalized_reaction',r/ref),**zcols('derivative',d)))
    write_csv(out/'field_reactions.csv',tab);write_csv(out/'reaction_density.csv',density)
    from .plots import reaction_plot,height_integrated_density_plot
    height_integrated_density_plot(density,out/'Fig3_from_fields.png')
    reaction_plot(result['RN']['La']/ref,result['RN']['Lb']/ref,out/'Fig7_from_fields.png','Computed from supplied field exports')
    write_json(out/'field_summary.json',{'provenance':'User-supplied field quadrature CSV; metadata declarations not independently verified','cancellation':cancellation(result['RN']['La']/ref,result['RN']['Lb']/ref),'window_symmetry_error':{k:relative_error(v['R'],v['L']) for k,v in result.items()}})
    return tab

FD_KEYS=['design','frequency_hz','window','delta','sign','s21_re','s21_im']
def fd(path,out,reaction_file=None):
    rs=rows(path,FD_KEYS); out=Path(out)
    for r in rs:
        check_frequency(r)
        if r['design'] not in DESIGNS or r['window'] not in ('L','R') or r['sign'] not in ('+','-'): raise ValueError('Invalid design/window/sign')
    ix=unique_map(rs,lambda r:(r['design'],r['window'],number(r['delta']),r['sign']))
    steps=(.01,.005,.0025);tab=[];checks=[];vals={}
    for k,w in itertools.product(DESIGNS,('L','R')):
        ds=[]
        for h in steps:
            z=central_difference(zrow(ix[k,w,h,'+'],'s21'),zrow(ix[k,w,h,'-'],'s21'),h);ds.append(z)
            tab.append(dict(design=k,window=w,delta=h,**zcols('derivative',z),magnitude=abs(z),phase_deg=angle(z)))
        vals[k,w]=ds[0]
        for i in range(2):
            checks.append(dict(design=k,window=w,coarse_step=steps[i],fine_step=steps[i+1],magnitude_difference_percent=100*abs(abs(ds[i])-abs(ds[i+1]))/abs(ds[i+1]),phase_difference_deg=phase_difference(ds[i],ds[i+1])))
    write_csv(out/'fd_derivatives.csv',tab);write_csv(out/'fd_step_halving.csv',checks)
    summary={'max_step_magnitude_difference_percent':max(x['magnitude_difference_percent'] for x in checks),'max_step_phase_difference_deg':max(x['phase_difference_deg'] for x in checks),'differential_common_percent':{k:100*differential_ratio(vals[k,'L'],vals[k,'R']) for k in DESIGNS}}
    if reaction_file:
        rec=unique_map(rows(reaction_file,['design','region','derivative_re','derivative_im']),lambda r:(r['design'],r['region']))
        ref=zrow(rec['ND','L'],'derivative'); cmp=[]
        for k,w in itertools.product(DESIGNS,('L','R')):
            a=zrow(rec[k,w],'derivative');b=vals[k,w]
            cmp.append(dict(design=k,window=w,**zcols('reaction_normalized',a/ref),**zcols('fd_normalized',b/ref),complex_relative_error_percent=100*relative_error(b,a)))
        write_csv(out/'reaction_fd_comparison.csv',cmp);summary['max_complex_relative_error_percent']=max(x['complex_relative_error_percent'] for x in cmp)
    write_json(out/'fd_summary.json',summary)
    return vals

FINITE_KEYS=['design','condition','frequency_hz','s21_re','s21_im','s11_re','s11_im','s22_re','s22_im','radiation_efficiency','eta_l_re','eta_l_im','eta_r_re','eta_r_im']
def finite(path,out,reaction_file=None):
    rs=rows(path,FINITE_KEYS);out=Path(out)
    for r in rs: check_frequency(r)
    ix=unique_map(rs,lambda r:(r['design'],r['condition']))
    p=read_json(ROOT/'inputs/paper.json');mat={k:material_eta(v['eps_real'],v['loss_tangent']) for k,v in p['materials'].items()}
    loads={'baseline':(0j,0j),**{q:tuple(mat[k] for k in v) for q,v in p['loads'].items()}}
    for k,q in itertools.product(DESIGNS,loads):
        r=ix[k,q]
        if not np.allclose([zrow(r,'eta_l'),zrow(r,'eta_r')],loads[q],rtol=0,atol=1e-9): raise ValueError(f'Wrong material loading for {k}/{q}')
    rec=unique_map(rows(reaction_file,['design','region','derivative_re','derivative_im']),lambda r:(r['design'],r['region'])) if reaction_file else None
    tab=[]
    for q in CONDITIONS:
        nd=abs(zrow(ix['ND',q],'s21')-zrow(ix['ND','baseline'],'s21'))
        if nd==0: raise ValueError('ND finite response is zero')
        for k in DESIGNS:
            r=ix[k,q];base=ix[k,'baseline'];z=zrow(r,'s21')-zrow(base,'s21');s11=zrow(r,'s11');s22=zrow(r,'s22');eff=number(r['radiation_efficiency'])
            line=dict(design=k,condition=q,**zcols('delta_s21',z),magnitude=abs(z),normalized_to_nd=abs(z)/nd,s11_db=float(db(s11)),s22_db=float(db(s22)),radiation_efficiency=eff,loaded_constraints_pass=bool(max(abs(s11),abs(s22))<=linear(-15) and .8<=eff<=1))
            if rec:
                h=remainder(zrow(r,'s21'),zrow(base,'s21'),zrow(rec[k,'L'],'derivative'),zrow(rec[k,'R'],'derivative'),*loads[q]);line.update(zcols('higher_order',h),higher_order_normalized=abs(h)/nd)
            tab.append(line)
    write_csv(out/'finite_responses.csv',tab);return tab

VNA_KEYS=['design','condition','cycle','frequency_hz','loaded_re','loaded_im','unloaded_re','unloaded_im']
def vna(path,out,expected_cycles=10,unloaded_bound=None):
    rs=rows(path,VNA_KEYS);out=Path(out);group={};cycle_rows=[]
    for r in rs:
        check_frequency(r);k=(r['design'],r['condition']);cy=number(r['cycle'])
        if k[0] not in DESIGNS or k[1] not in CONDITIONS or cy!=int(cy) or cy<1: raise ValueError('Invalid design/condition/cycle')
        z=zrow(r,'loaded')-zrow(r,'unloaded'); group.setdefault(k,{})
        if cy in group[k]: raise ValueError(f'Duplicate cycle {k}/{cy}')
        group[k][cy]=abs(z);cycle_rows.append(dict(design=k[0],condition=k[1],cycle=int(cy),**zcols('delta_s21',z),magnitude=abs(z)))
    tab=[]
    for q in CONDITIONS:
        nd=np.mean(list(group['ND',q].values()))
        if nd<=0: raise ValueError('Zero ND mean')
        for k in DESIGNS:
            a=np.array(list(group[k,q].values()));n=len(a)
            if n!=expected_cycles or n<2 or set(group[k,q])!=set(range(1,expected_cycles+1)): raise ValueError(f'{k}/{q}: need cycles 1..{expected_cycles}')
            mean=float(a.mean());sd=float(a.std(ddof=1));line=dict(design=k,condition=q,n=n,mean=mean,sample_sd=sd,normalized_mean=mean/nd,normalized_sd=sd/nd,suppression=nd/mean if mean else '')
            if unloaded_bound is not None:
                if unloaded_bound<=0: raise ValueError('Positive unloaded bound required')
                line['mean_to_unloaded_bound']=mean/unloaded_bound
            tab.append(line)
    write_csv(out/'vna_cycle_magnitudes.csv',cycle_rows);write_csv(out/'vna_summary.csv',tab)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,4));
    for i,k in enumerate(DESIGNS):
        a=[r for r in tab if r['design']==k];ax.errorbar(np.arange(2)+(i-1)*.15,[r['normalized_mean'] for r in a],yerr=[r['normalized_sd'] for r in a],marker=['s','^','o'][i],label=k,linestyle='none',capsize=4)
    ax.set_xticks([0,1],CONDITIONS);ax.set_ylabel('Normalized mean ± 1 sample SD');ax.set_title('Supplied VNA cycles');ax.legend();ax.grid(alpha=.2);save(fig,out/'Fig10_vna_from_raw.png');return tab
