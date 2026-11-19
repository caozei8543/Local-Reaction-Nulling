"""Synthetic unit-test fixtures ONLY, in temporary directories, never research data."""
import unittest,tempfile,itertools,subprocess,sys
from pathlib import Path
import numpy as np
from rn.core import *
from rn.io import *
from rn import raw,solver,grid,touchstone,reproduce,sweeps

class TestRN(unittest.TestCase):
    def setUp(self): self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def test_bilinear_not_hermitian(self):
        a=np.array([[1+1j,2j,3]]);b=np.array([[2j,1j,2]])
        r=reaction(a,b,[2]);self.assertAlmostEqual(r,4+4j);self.assertNotEqual(r,np.vdot(a,b)*2)
    def test_normalization(self):
        r=2+3j;d=derivative_from_reaction(r)
        self.assertAlmostEqual(derivative_from_reaction(2*r,phasor='peak'),d)
        self.assertAlmostEqual(derivative_from_reaction(6*r,p1_w=4,p2_w=9),d)
        self.assertAlmostEqual(derivative_from_reaction(r),-1j*np.pi*F0*EPS0*r)
    def test_central_complex_and_symmetry(self):
        d=.12+.03j;s0=.02+.01j;h=.01
        f=lambda t:s0+d*t+.3j*t*t
        self.assertAlmostEqual(central_difference(f(h),f(-h),h),d)
        self.assertAlmostEqual(differential_ratio(d,d),0)
        self.assertAlmostEqual(delta_linear(d,d,.1,-.1),0)
        self.assertAlmostEqual(worst_case_first_order(d,d,.1),.2*abs(d))
    def test_paper(self):
        s=reproduce.run(self.p/'reported');self.assertAlmostEqual(s['cancellation']['sum_magnitude'],.02003916415422554)
        self.assertAlmostEqual(s['cancellation']['cres'],.01391603918553)
        self.assertAlmostEqual(s['max_fd_relative_error_percent'],5.03902896775)
        self.assertAlmostEqual(s['max_differential_common_percent'],1.07128827181)
        self.assertAlmostEqual(s['fd_suppression']['left'],47.38095238095)
        self.assertAlmostEqual(s['fd_suppression']['right'],47.93269230769)
        self.assertAlmostEqual(s['measured_suppression']['common'],20.83333333333)
        self.assertAlmostEqual(s['measured_suppression']['unequal'],18.18181818182)
        self.assertAlmostEqual(s['rn_total_efficiency_percent'],85.76995120757)
    def test_touchstone(self):
        for typ,values in [('RI','0.1 0 0.01 0.02 0.01 0.02 0.1 0'),('MA','0.1 0 0.05 90 0.05 90 0.1 0'),('DB','-20 0 -40 90 -40 90 -20 0')]:
            f=self.p/f'{typ}.s2p';f.write_text(f'# GHz S {typ} R 50\n5.8 {values}\n');expected={'RI':.01+.02j,'MA':.05j,'DB':.01j}[typ];self.assertAlmostEqual(touchstone.s21_at_f0(f),expected)
    def test_vna(self):
        lines=[]
        for k,q,i in itertools.product(DESIGNS,CONDITIONS,range(1,11)):
            scale={'ND':1,'DLO':.5,'RN':.1}[k];z=.1+1j*.01*i*scale
            lines.append(dict(design=k,condition=q,cycle=i,frequency_hz=F0,**zcols('loaded',z),**zcols('unloaded',.1)))
        f=self.p/'v.csv';write_csv(f,lines);v=raw.vna(f,self.p/'v');r=next(x for x in v if x['design']=='RN');self.assertAlmostEqual(r['normalized_mean'],.1);self.assertAlmostEqual(r['sample_sd'],np.std(np.arange(1,11)*.001,ddof=1))
        write_csv(f,lines[:-1]);self.assertRaises(ValueError,raw.vna,f,self.p/'bad')
    def test_optimizer_budget_and_constraints(self):
        p=np.mean(BOUNDS,axis=1)
        def toy(x,l,r,need):
            # Analytic test function, not an electromagnetic model.
            return dict(s11=.05,s22=.05,s21=.05+.01*(x[0]-1.2)**2+(.02+.01j)*(l+r),radiation_efficiency=.9,reaction_l=complex(x[0]-1.2,.01),reaction_r=complex(x[0]-1.2,.01))
        for mode in DESIGNS:
            r=solver.optimize_one(mode,toy,p,1+0j);self.assertLessEqual(r['evaluations'],30);self.assertLessEqual(r['iterations'],25);self.assertTrue(all(x['feasible'] for x in r['history']))
    def test_grid_and_merge(self):
        m={'regions_mm':{k:{'min':[i*10,0,0],'max':[i*10+9,3,2],'cells':[3,2,2]} for i,k in enumerate(raw.REGIONS)}};f=self.p/'g.json';write_json(f,m);g=self.p/'g.csv';a=grid.generate(f,g);self.assertEqual(len(a),48);self.assertAlmostEqual(sum(r['weight_m3'] for r in a),216e-9)
        data=[dict(x_m=r['x_m'],y_m=r['y_m'],z_m=r['z_m'],**{f'E{ax}_{c}':1 if c=='re' else 0 for ax in 'xyz' for c in ('re','im')}) for r in a];q=self.p/'p.csv';write_csv(q,data);self.assertEqual(len(grid.merge(g,q,q,'ND',self.p/'m.csv')),48)
    def test_pipeline(self):
        # Entire end-to-end pipeline, synthetic algebraic inputs used ONLY for unit testing.
        p=self.p;f=[];fd=[];finite=[];vna=[];sw=[];ds={};m=read_json(ROOT/'inputs/templates/field_metadata.json');m['phasor']='rms';write_json(p/'meta.json',m)
        for i,k in enumerate(DESIGNS):
            scale=(1,.5,.1)[i];rvals=[1+.2j,(-1+.1*scale)-.2j]*2;dl=derivative_from_reaction(sum(rvals[:2]));ds[k]=dl;s0=.05+.001j
            for j,q in enumerate(raw.REGIONS):
                r=dict(design=k,frequency_hz=F0,region=q,x_m=j*.01,y_m=0,z_m=.001,weight_m3=54e-9,g=1)
                for ax in 'xyz':r.update(**zcols(f'E1{ax}',1 if ax=='x' else 0),**zcols(f'E2{ax}',rvals[j]/54e-9 if ax=='x' else 0))
                f.append(r)
            for w,h,sign in itertools.product(('L','R'),(.01,.005,.0025),('+','-')):
                z=s0+dl*h*(1 if sign=='+' else -1);fd.append(dict(design=k,frequency_hz=F0,window=w,delta=h,sign=sign,**zcols('s21',z)))
            for q,l,r in [('baseline',0j,0j),('common',material_eta(1.06,.0008),material_eta(1.06,.0008)),('unequal',material_eta(1.06,.0008),material_eta(1.106,.0016))]:
                z=s0+dl*(l+r);finite.append(dict(design=k,condition=q,frequency_hz=F0,**zcols('s21',z),**zcols('s11',.05),**zcols('s22',.05),radiation_efficiency=.9,**zcols('eta_l',l),**zcols('eta_r',r)))
            for q,j in itertools.product(CONDITIONS,range(1,11)):
                vna.append(dict(design=k,condition=q,cycle=j,frequency_hz=F0,**zcols('loaded',s0+1j*.001*j*scale),**zcols('unloaded',s0)))
            for et,fr in itertools.product((-.1,0,.1),(5.5e9,5.8e9,6.1e9)):
                sw.append(dict(design=k,source='simulation',eta_l=et,eta_r=0,frequency_hz=fr,**zcols('s11',.1),**zcols('s21',.05+((fr-F0)/1e9-et)**2)))
        for key,data in [('fields',f),('fd',fd),('finite',finite),('vna',vna),('sweeps',sw)]:write_csv(p/(key+'.csv'),data)
        write_json(p/'manifest.json',{'fields_csv':'fields.csv','field_metadata_json':'meta.json','fd_csv':'fd.csv','finite_csv':'finite.csv','vna_csv':'vna.csv','sweeps_csv':'sweeps.csv'})
        r=subprocess.run([sys.executable,'-m','rn','pipeline','--manifest',str(p/'manifest.json'),'--out',str(p/'out')],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr);self.assertTrue((p/'out'/'Fig10_from_raw.png').exists());self.assertLess(read_json(p/'out'/'fd'/'fd_summary.json')['max_complex_relative_error_percent'],1e-7)
    def test_missing_data_not_faked(self):
        self.assertRaises(ValueError,raw.vna,ROOT/'inputs/templates/vna.csv',self.p/'v')

if __name__=='__main__':unittest.main()
