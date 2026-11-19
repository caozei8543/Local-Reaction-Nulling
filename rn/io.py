"""Strict CSV/JSON input schemas. Missing source data are errors, not generated replacements."""
import csv,json
from pathlib import Path
import numpy as np
from .core import F0
ROOT=Path(__file__).resolve().parents[1]

def read_json(path):
    with open(path,encoding='utf-8') as f: return json.load(f)

def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8') as f: json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False); f.write('\n')

def rows(path, required):
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or not set(required)<=set(reader.fieldnames): raise ValueError(f'{path}: missing columns {set(required)-set(reader.fieldnames or [])}')
        result=list(reader)
    if not result: raise ValueError(f'{path}: no data; supply genuine exports')
    return result

def write_csv(path,data,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    data=list(data)
    fields=fields or (list(data[0]) if data else None)
    if fields is None: raise ValueError('No columns specified')
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)

def number(value):
    v=float(value)
    if not np.isfinite(v): raise ValueError(f'Nonfinite number: {value}')
    return v

def zrow(row,name): return complex(number(row[name+'_re']),number(row[name+'_im']))

def zcols(name,z): return {name+'_re':float(np.real(z)),name+'_im':float(np.imag(z))}

def check_frequency(row):
    if abs(number(row['frequency_hz'])-F0)>1: raise ValueError('Expected exactly 5.8 GHz (1 Hz tolerance); no retuning/interpolation')

def unique_map(records,key):
    result={}
    for r in records:
        k=key(r)
        if k in result: raise ValueError(f'Duplicate input: {k}')
        result[k]=r
    return result
