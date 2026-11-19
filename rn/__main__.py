import argparse,json
from pathlib import Path
from . import reproduce,raw,sweeps,solver,touchstone,grid
from .io import ROOT

def main():
    parser=argparse.ArgumentParser(description='Local reaction nulling: source-grounded reproduction and real-data analysis')
    sub=parser.add_subparsers(dest='cmd',required=True)
    def cmd(name,help):
        p=sub.add_parser(name,help=help);p.add_argument('--out',type=Path,required=True);return p
    p=cmd('reproduce','Reconstruct paper tables and plots');p.add_argument('--paper',type=Path,default=ROOT/'inputs/paper.json')
    p=cmd('fields','Integrate true nominal E-field exports');p.add_argument('--csv',required=True);p.add_argument('--metadata',required=True)
    for n in ('fd','finite'):
        p=cmd(n,'Process true complex simulation data');p.add_argument('--csv',required=True);p.add_argument('--reactions')
    p=cmd('vna','Process ten paired complex VNA cycles');p.add_argument('--csv',required=True);p.add_argument('--unloaded-bound',type=float)
    p=cmd('sweeps','Analyze frequency sweeps');p.add_argument('--csv',required=True)
    p=cmd('touchstone','Convert loaded/unloaded .s2p manifest into VNA cycle CSV');p.add_argument('--manifest',required=True)
    p=cmd('optimize','Run reference constrained optimizer through external HFSS worker');p.add_argument('--solver',required=True);p.add_argument('--initial',required=True)
    p=cmd('collect','Run fixed-geometry FD and finite loads through external HFSS worker');p.add_argument('--solver',required=True);p.add_argument('--geometry')
    p=cmd('grid','Generate common HFSS quadrature grid from actual CAD bounds');p.add_argument('--config',required=True)
    p=cmd('merge-fields','Merge common grid and two nominal field exports');p.add_argument('--grid',required=True);p.add_argument('--port1',required=True);p.add_argument('--port2',required=True);p.add_argument('--design',choices=['ND','DLO','RN'],required=True)
    p=cmd('pipeline','Run complete raw CSV analysis from one dataset manifest');p.add_argument('--manifest',required=True)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if a.cmd=='reproduce': print(json.dumps(reproduce.run(a.out,a.paper),indent=2,ensure_ascii=False))
    elif a.cmd=='fields':raw.fields(a.csv,a.metadata,a.out)
    elif a.cmd=='fd':raw.fd(a.csv,a.out,a.reactions)
    elif a.cmd=='finite':raw.finite(a.csv,a.out,a.reactions)
    elif a.cmd=='vna':raw.vna(a.csv,a.out,unloaded_bound=a.unloaded_bound)
    elif a.cmd=='sweeps':sweeps.sweeps(a.csv,a.out)
    elif a.cmd=='touchstone':touchstone.convert_manifest(a.manifest,a.out/'vna_cycles.csv')
    elif a.cmd=='optimize':solver.optimize_all(a.solver,a.initial,a.out)
    elif a.cmd=='collect':solver.collect(a.solver,a.out,a.geometry)
    elif a.cmd=='grid':grid.generate(a.config,a.out/'quadrature_grid.csv')
    elif a.cmd=='merge-fields':grid.merge(a.grid,a.port1,a.port2,a.design,a.out/'fields.csv')
    elif a.cmd=='pipeline':
        path=Path(a.manifest);m=json.loads(path.read_text(encoding='utf-8'));get=lambda k:path.parent/m[k]
        raw.fields(get('fields_csv'),get('field_metadata_json'),a.out/'fields');r=a.out/'fields'/'field_reactions.csv'
        raw.fd(get('fd_csv'),a.out/'fd',r);raw.finite(get('finite_csv'),a.out/'finite',r);raw.vna(get('vna_csv'),a.out/'vna',unloaded_bound=m.get('unloaded_bound'));sweeps.sweeps(get('sweeps_csv'),a.out/'sweeps')
        from .plots import comparison_plot
        from .io import rows,zrow
        from .core import DESIGNS,CONDITIONS
        fdrows=rows(a.out/'fd'/'reaction_fd_comparison.csv',['design','window','fd_normalized_re','fd_normalized_im']);vr=rows(a.out/'vna'/'vna_summary.csv',['design','condition','normalized_mean','normalized_sd'])
        data={};errors={}
        for k in DESIGNS:
            data[k]={f'fd_{w}_mag_deg':[abs(zrow(next(r for r in fdrows if r['design']==k and r['window']==w[0].upper()),'fd_normalized')),0] for w in ('left','right')}
            data[k]['meas_finite_ratios']={q:float(next(r for r in vr if r['design']==k and r['condition']==q)['normalized_mean']) for q in CONDITIONS}
            errors[k]=[float(next(r for r in vr if r['design']==k and r['condition']==q)['normalized_sd']) for q in CONDITIONS]
        comparison_plot(data,a.out/'Fig10_from_raw.png',errors,'Supplied FD and VNA exports')
    print(f'Completed {a.cmd}: {a.out.resolve()}')
if __name__=='__main__':main()
