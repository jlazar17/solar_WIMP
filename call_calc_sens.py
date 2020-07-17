import numpy as np

from calc_sens import SensCalc
from controls import mcfiles, chs, ms

def initialize_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--syst_set')
    args=parser.parse_args()
    return args

def main(sigflux, syst_set):
    sc = SensCalc(sigflux)
    for mcpath in mcfiles:
        if syst_set in mcpath:
            sc.add_mus(mcpath)
            sc.add_mub(mcpath)
    sc.set_sens()
    return sc.sens


if __name__=='__main__':
    args = initialize_args()
    sens_dict = {}
    flux_types = [f'ch{ch}-m{m}' for ch in chs
                                 for m in ms]
    for ft in flux_types:
        try:
            sens_dict[ft] = main(ft, args.syst_set)
        except:
            print(ft)
    print(sens_dict)
    print('/data/user/jlazar/solar_WIMP/data/sens_dicts/%s_sens_dict.npy' % args.syst_set)
    np.save('/data/user/jlazar/solar_WIMP/data/sens_dicts/%s_sens_dict.npy' % args.syst_set, sens_dict, allow_pickle=True)
