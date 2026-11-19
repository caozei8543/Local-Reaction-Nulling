"""Exact response lookup for an archive of actual HFSS exports; never simulates.
Usage via solver command: python tools/export_archive_worker.py real_index.json request.json response.json
"""
from pathlib import Path
import json,sys

def main():
    if len(sys.argv)!=4:raise SystemExit('usage: export_archive_worker.py INDEX REQUEST RESPONSE')
    index,req,out=map(Path,sys.argv[1:]);r=json.loads(req.read_text(encoding='utf-8'));items=json.loads(index.read_text(encoding='utf-8'))
    match=[x for x in items if x.get('request')==r and x.get('status')=='ok' and x.get('provenance')=='full_wave_export']
    if len(match)!=1:raise SystemExit('No unique exact real HFSS export for request; run actual solver; refusing to fabricate results.')
    out.write_text(json.dumps(match[0],indent=2,allow_nan=False),encoding='utf-8')
if __name__=='__main__':main()
