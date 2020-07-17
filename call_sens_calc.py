import numpy as np


from calc_sens import SenseCalc
from controls import mcfiles, ms, chs

def initialize_parser():
    import argparse
    parser = argpare.ArgumentParser()
    parser.add_argument('--syst_set', type=str)
    return parser.parse_args()

def main(sigflux, syst_set)
    sc = SensCalc(sigflux)
    for mcpath in mcfiles:
        if syst_set in mcpath:
            sc.add_mus(macpath)
            sc.add_mub(macpath)
    sc.set_sens
    return sc.sens

if __name__=='__main__':
    args = initialize_parser()
    sigfluxes = ['ch%d-m%d' % (ch,m) for ch in chs
                                     for m in ms]
    d = {}
    for sf in sigfluxes:
        d[sf] = main(sf, args.syst_set)

    print(d)
