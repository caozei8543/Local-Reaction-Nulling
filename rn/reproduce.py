"""Reconstruct only metrics mathematically determined by the two PDFs."""
from pathlib import Path
import numpy as np
from .core import *
from .io import *

def run(out,source=ROOT/'inputs/paper.json'):
    p=read_json(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ref=complex(*p['reference_derivative_re_im']); a,b=[complex(*x) for x in p['rn_subregion_reactions_normalized_re_im']]
    c=cancellation(a,b); mat={k:material_eta(v['eps_real'],v['loss_tangent']) for k,v in p['materials'].items()}
    loads={k:tuple(mat[v] for v in vals) for k,vals in p['loads'].items()}
    geom=[];nom=[];sens=[];absolute=[];check=[];lin=[];finite=[];means=[]
    d=p['designs']
    for i,n in enumerate(p['parameter_order']):
        geom.append({'parameter':n,'lower_mm':p['bounds_mm'][i][0],'upper_mm':p['bounds_mm'][i][1],**{k:d[k]['geometry_mm'][i] for k in DESIGNS}})
    for k in DESIGNS:
        v=d[k];r,l,h=[polar(*v[n]) for n in ('reaction_mag_deg','fd_left_mag_deg','fd_right_mag_deg')]
        nom.append(dict(design=k,s21_db=v['s21_db'],s11_db=v['s11_db'],radiation_efficiency_percent=100*v['radiation_efficiency'],measured_isolation_db=v['meas_isolation_db']))
        sens.append(dict(design=k,reaction_magnitude=abs(r),reaction_phase_deg=angle(r),fd_left_magnitude=abs(l),fd_left_phase_deg=angle(l),fd_right_magnitude=abs(h),fd_right_phase_deg=angle(h),**{'sim_'+q:v['sim_finite_ratios'][q] for q in CONDITIONS}))
        absolute.append(dict(design=k,**zcols('reaction',r*ref),**zcols('fd_left',l*ref),**zcols('fd_right',h*ref)))
        check.append(dict(design=k,error_left_percent=100*relative_error(l,r),error_right_percent=100*relative_error(h,r),differential_common_percent=100*differential_ratio(l,h)))
        lin.append(dict(design=k,**{q:abs(delta_linear(r*ref,r*ref,*loads[q])) for q in CONDITIONS}))
        finite.append(dict(design=k,**{'sim_'+q:v['sim_finite_ratios'][q] for q in CONDITIONS},**{'meas_'+q:v['meas_finite_ratios'][q] for q in CONDITIONS}))
        means.append(dict(design=k,**{q:v['meas_finite_ratios'][q]*p['measurement']['nd_means'][q] for q in CONDITIONS}))
    for name,tab in [('Table_I_S1',geom),('Table_II_S2',nom),('Table_III_S3',sens),('Table_S4',absolute),('Table_S5',check),('Table_S6',lin),('Table_S7',finite),('Table_S8',means)]: write_csv(out/(name+'.csv'),tab)
    rn=d['RN'];nd=d['ND']; efficiency=100*total_efficiency(rn['radiation_efficiency'],linear(rn['s11_db']),linear(rn['s21_db']))
    summary={'provenance':p['provenance'],'cancellation':c,'material_eta_re_im':{k:pair(v) for k,v in mat.items()},'reference_reaction_SI_re_im':pair(ref/derivative_from_reaction(1)), 'fd_suppression':{q:nd['fd_'+q+'_mag_deg'][0]/rn['fd_'+q+'_mag_deg'][0] for q in ('left','right')},'measured_suppression':{q:1/rn['meas_finite_ratios'][q] for q in CONDITIONS},'simulated_finite_suppression':{q:1/rn['sim_finite_ratios'][q] for q in CONDITIONS},'rn_total_efficiency_percent':efficiency,'rn_fractional_bandwidth_percent':100*p['rn_nominal_bandwidth_mhz']*1e6/F0,'rn_mean_to_reported_unloaded_upper_bound':{q:means[-1][q]/p['measurement']['unloaded_variation_upper_bound'] for q in CONDITIONS},'max_fd_relative_error_percent':max(max(x['error_left_percent'],x['error_right_percent']) for x in check),'max_differential_common_percent':max(x['differential_common_percent'] for x in check),'notes':['Table S4/S5/S6/S8 are algebraic reconstructions from rounded reported values. Last digits can differ.','FD ratios use actual rounded ND FD numerator; finite-load ratios use their own ND denominator.','Fig. 10 reconstructed plot has no error bars: individual cycles and SD are not in the PDFs.','No numerical frequency trajectories, field samples, optimization history, or VNA cycles were reconstructed.']}
    write_json(out/'summary.json',summary)
    from .plots import reaction_plot,comparison_plot
    reaction_plot(a,b,out/'Fig7_S1_reconstruction.png')
    comparison_plot(d,out/'Fig10_reported_ratios.png')
    return summary
