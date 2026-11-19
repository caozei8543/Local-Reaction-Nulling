"""Solver boundary: real full-wave solutions are provided by the user's HFSS worker.
No electromagnetic data are synthesized here. A deterministic local search is supplied,
not the unknown optimizer or history used in the manuscript.
"""
import json,subprocess,tempfile,itertools
from pathlib import Path
import numpy as np
from .core import *
from .io import *

class ExternalSolver:
    def __init__(self,config):
        self.config=read_json(config);self.command=self.config['command'];self.timeout=float(self.config.get('timeout_seconds',3600));self.calls=0
        if not isinstance(self.command,list) or not self.command or not all(isinstance(x,str) for x in self.command): raise ValueError('command must be an argv list, not shell text')
    def __call__(self,p,eta_l=0j,eta_r=0j,need_reaction=False):
        req={'geometry_mm':dict(zip(PARAMETERS,map(float,p))),'eta_l_re_im':pair(eta_l),'eta_r_re_im':pair(eta_r),'frequency_hz':F0,'reference_ohm':50,'need_reaction':need_reaction}
        with tempfile.TemporaryDirectory() as t:
            a=Path(t)/'request.json';b=Path(t)/'response.json';write_json(a,req)
            self.calls+=1;subprocess.run([*self.command,str(a),str(b)],check=True,timeout=self.timeout)
            r=read_json(b)
        if r.get('status')!='ok' or r.get('request')!=req: raise ValueError('Worker must acknowledge exact request')
        if r.get('provenance')!='full_wave_export': raise ValueError('Expected actual full-wave export provenance')
        z={k:complex(*r[k+'_re_im']) for k in ('s11','s21','s22')};z['radiation_efficiency']=number(r['radiation_efficiency'])
        if need_reaction:
            # Worker supplies unit-incident-power RMS normalized field integrals in SI.
            z['reaction_l']=complex(*r['reaction_l_SI_re_im']);z['reaction_r']=complex(*r['reaction_r_SI_re_im'])
        if not all(np.isfinite(v) for v in z.values()) or not 0<=z['radiation_efficiency']<=1: raise ValueError('Invalid solution response')
        return z

def optimize_one(mode,solve,initial,reference=None,max_iterations=25,max_evaluations=30):
    """Feasible-only coordinate poll, same initial point/bounds/budget for all objectives.
    This is a reference implementation, not a reconstruction of the unreported algorithm.
    """
    if mode not in DESIGNS or max_iterations<0 or max_evaluations<1: raise ValueError('Invalid optimization settings')
    if mode=='RN' and (reference is None or abs(reference)==0): raise ValueError('RN requires fixed optimized ND reference reaction')
    p=np.asarray(initial,float); history=[]; cache={}; solver_calls=0
    if p.shape!=(7,) or not np.isfinite(p).all() or np.any(p<BOUNDS[:,0]) or np.any(p>BOUNDS[:,1]): raise ValueError('Initial geometry out of bounds')
    def call(x,l=0j,r=0j,need=False):
        nonlocal solver_calls
        solver_calls+=1;return solve(x,l,r,need)
    def evaluate(x):
        key=tuple(np.round(x,12))
        if key in cache:return cache[key]
        if len(history)>=max_evaluations:return (False,float('inf'))
        n=call(x,need=mode=='RN');ok=feasible(x,n['s11'],n['s21'],n['s22'],n['radiation_efficiency']);obj=None
        if ok:
            if mode=='ND':obj=abs(n['s21'])**2
            elif mode=='RN':obj=abs(n['reaction_l']/reference)**2
            else:
                obj=sum(abs(call(x,t if w=='L' else 0,t if w=='R' else 0)['s21']-n['s21'])**2 for w,t in itertools.product(('L','R'),(-.1,-.05,.05,.1)))
        history.append(dict(evaluation=len(history)+1,geometry_mm=list(map(float,x)),feasible=ok,objective=obj));cache[key]=(ok,obj if ok else float('inf'));return cache[key]
    valid,best=evaluate(p)
    if not valid: raise ValueError('Common initial geometry must satisfy the nominal feasible set; supply the real Stage I solution')
    widths=(BOUNDS[:,1]-BOUNDS[:,0])*.1;iterations=0
    for it in range(max_iterations):
        if len(history)>=max_evaluations:break
        iterations=it+1; improved=False
        for i,sign in itertools.product(range(7),(1,-1)):
            q=p.copy();q[i]=np.clip(q[i]+sign*widths[i],*BOUNDS[i]);ok,v=evaluate(q)
            if ok and v<best:
                p,best=q,v;improved=True;break
            if len(history)>=max_evaluations:break
        if not improved:widths*=.5
        if widths.max()<1e-7:break
    return {'mode':mode,'geometry_mm':dict(zip(PARAMETERS,map(float,p))),'objective':best,'iterations':iterations,'evaluations':len(history),'worker_requests':solver_calls,'history':history,'algorithm':'Reference coordinate poll; not manuscript optimizer history'}

def optimize_all(config,initial_file,out):
    out=Path(out);s=ExternalSolver(config);initial=read_json(initial_file)['initial_geometry_mm'];p=[initial[k] for k in PARAMETERS]
    nd=optimize_one('ND',s,p);ref=s([nd['geometry_mm'][k] for k in PARAMETERS],0j,0j,True)['reaction_l'];dlo=optimize_one('DLO',s,p);rn=optimize_one('RN',s,p,ref)
    result={'ND':nd,'DLO':dlo,'RN':rn,'nd_reference_reaction_SI_re_im':pair(ref),'total_worker_requests_including_reference':s.calls};write_json(out/'optimization.json',result);return result

def collect(config,out,geometry_file=None):
    """Run real worker to collect FD, finite-load and reaction data for fixed designs.
    Defaults to Table I; geometry override uses optimize_all result. No retuning.
    """
    from .raw import fd,finite
    out=Path(out);p=read_json(ROOT/'inputs/paper.json');s=ExternalSolver(config);g=read_json(geometry_file) if geometry_file else None
    mat={k:material_eta(v['eps_real'],v['loss_tangent']) for k,v in p['materials'].items()};load={'baseline':(0j,0j),**{q:tuple(mat[v] for v in a) for q,a in p['loads'].items()}};fdr=[];fr=[];rr=[]
    for k in DESIGNS:
        x=[g[k]['geometry_mm'][n] for n in PARAMETERS] if g else p['designs'][k]['geometry_mm']
        baseline=s(x,0j,0j,True)
        for w in ('L','R'):
            r=baseline['reaction_'+w.lower()];rr.append(dict(design=k,region=w,**zcols('reaction',r),**zcols('derivative',derivative_from_reaction(r))))
        for w,h,sign in itertools.product(('L','R'),p['fd_steps'],('+','-')):
            t=h if sign=='+' else -h;v=s(x,t if w=='L' else 0,t if w=='R' else 0)
            fdr.append(dict(design=k,window=w,delta=h,sign=sign,frequency_hz=F0,**zcols('s21',v['s21'])))
        for q,(l,r) in load.items():
            v=baseline if q=='baseline' else s(x,l,r)
            fr.append(dict(design=k,condition=q,frequency_hz=F0,**zcols('s21',v['s21']),**zcols('s11',v['s11']),**zcols('s22',v['s22']),radiation_efficiency=v['radiation_efficiency'],**zcols('eta_l',l),**zcols('eta_r',r)))
    write_csv(out/'solver_reactions.csv',rr);write_csv(out/'solver_fd.csv',fdr);write_csv(out/'solver_finite.csv',fr)
    fd(out/'solver_fd.csv',out,out/'solver_reactions.csv');finite(out/'solver_finite.csv',out,out/'solver_reactions.csv');write_json(out/'solver_collection.json',{'worker_requests':s.calls,'geometry_source':'supplied optimizer output' if g else 'paper Table I','provenance':'external full-wave worker; not independently certified'})
