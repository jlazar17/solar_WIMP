import numpy as np
import argparse
from os import path

from path_gen import PathGen
from controls_2 import mc_list, mc_types, default_ms, default_chs, sens_dict
from calc_signal_flux import FluxCalculator
from calc_mc_flux import create_mc_fluxmaker
from calc_e_d_theta_hist_submit import make_dag_file
from calc_sensitivity import SensitivityCalculator

def initialize_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcoptions',
                        type=str,
                        default='1111111111',
                        help='Which systematic sets to process.'
                       )
    parser.add_argument('--masses',
                        nargs='+',
                        default=default_ms,
                        help='Which masses to include.'
                       )
    parser.add_argument('--channels',
                        nargs='+',
                        default=default_chs,
                       )
    parser.add_argument('--recompute',
                        default=False,
                        action='store_true'
                       )
    parser.add_argument('--edtheta_opts',
                        default='00',
                        type=str
                       )

    args = parser.parse_args()
    return args

class Analysis():

    def __init__(self,args ):
        self._mcoptions = args.mcoptions
        self.channels   = args.channels
        self.masses     = args.masses
        self.recompute  = args.recompute
        self.edtheta_opts = args.edtheta_opts
    
        self.mcfile_list = self.make_mcfile_list()
        self.fluxtype_list = self.make_flux_list()
        self.syst_type_list   = [mc_types[i] for i in range(len(mc_types)) if int(self._mcoptions[i])]
        print(self.syst_type_list)

    def make_mcfile_list(self):
        mcfile_list = []
        for ls, opt in zip(mc_list, self._mcoptions):
            if int(opt):
                while len(ls)!=0:
                    mcfile_list.append(ls.pop())
        return mcfile_list

    def make_flux_list(self):
        fluxtype_list = ['ch%d-m%d' % (ch,m) for ch in self.channels
                                              for m in self.masses]
        fluxtype_list.append('conv-numu')
        fluxtype_list.append('solar-atm')
        return fluxtype_list

    def compute_flux_at_earth(self):
        for ch in self.channels:
            for m in self.masses:
                flux_path = "/data/user/jlazar/solar_WIMP/data/charon_fluxes/ch%d-m%d_1AU_BRW_dn_dz.npy" % (ch, m)
                print(flux_path)
                if not path.exists(flux_path) or self.recompute:
                    print('computing flux for channel %d, mass %d...'% (ch, m))
                    fc = FluxCalculator(ch, m, '1AU', 'BRW')
                    fc.save_flux()

    def compute_mc_fluxes(self):
        for mcfile in self.mcfile_list:
            pg = PathGen(mcfile)
            for fluxtype in self.fluxtype_list:
                if not path.exists(pg.get_mc_dn_dz_path(fluxtype)) or self.recompute:
                    print(pg.get_mc_dn_dz_path(fluxtype))
                    fluxmaker = create_mc_fluxmaker(mcfile, fluxtype)
                    fluxmaker.initialize_nuSQuIDS()
                    fluxmaker.save_flux()

    def make_e_dtheta_list(self):
        e_dtheta_list = []
        for mcfile in self.mcfile_list:
            pg = PathGen(mcfile)
            for fluxtype in self.fluxtype_list:
                if not path.exists("%s_%s.npy" % (pg.get_e_d_theta_path(fluxtype),self.edtheta_opts)) or self.recompute:
                    e_dtheta_list.append((mcfile, fluxtype, self.edtheta_opts))
        return e_dtheta_list

    def make_edtheta_dag(self):
        missing_hists = self.make_e_dtheta_list()
        print(len(missing_hists))
        if len(missing_hists)==0:
            print('All histograms made.')
        else:
            make_dag_file(missing_hists)
        pass

    def calc_sensitivity(self):
        sensitivities = np.full((len(self.channels), len(self.masses), len(self.syst_type_list)), np.nan)
        for ic, ch in enumerate(self.channels):
            for im, m in enumerate(self.masses):
                for ist, syst_type in enumerate(self.syst_type_list):
                    try:
                        sc = SensitivityCalculator(ch, m, syst_type, self.edtheta_opts)
                        sensitivities[ic][im][ist] = sc.calc_sens()
                    except IOError:
                        print('Unable to load necessary histograms for ch %d, m %d, %s, %s' % (ch, m, syst_type, self.edtheta_opts))
        print(sensitivities)


if __name__=='__main__':
    args = initialize_args()
    analysis = Analysis(args)
    analysis.compute_flux_at_earth()
    analysis.compute_mc_fluxes()
    #analysis.make_edtheta_dag()
    analysis.calc_sensitivity()
