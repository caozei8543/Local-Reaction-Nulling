"""Convenient entrypoint: python run_reproduction.py"""
from rn.reproduce import run
from rn.io import ROOT
if __name__=='__main__':
    run(ROOT/'outputs/reported')
    print('Reproduced paper-summary outputs in outputs/reported/')
