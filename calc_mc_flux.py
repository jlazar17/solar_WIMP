import time
import h5py
import numpy as np

import nuSQUIDSpy as nsq
from path_gen import PathGen
from physicsconstants import PhysicsConstants
pc = PhysicsConstants()

def initialize_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcfile',
                        type=str,
                        help='path to mcfile to be used'
                       )
    parser.add_argument('--fluxtype',
                        type=str,
                        help='fluxtype to be used, e.g. conv-numu or ch5-m1000'
                       )
    parser.add_argument('--options',
                        type=str,
                        default='00'
                       )
    
    args = parser.parse_args()
    return args

class BaseMCFluxMaker:

    def __init__(self, mcpath, fluxtype):
        self.mcpath   = mcpath
        self.fluxtype = fluxtype
        self.mcdata   = self.read_mc()

    def read_mc(self):
        mcf = h5py.File(self.mcpath, 'r')
        keys2grab = ['NuEnergy', 'NuZenith', 'PrimaryType']
        mcdata = np.array(zip(mcf['NuEnergy'][()]['value'],
                              np.cos(mcf['NuZenith'][()]['value']),
                              mcf['PrimaryType'][()]['value']
                             ),
                          dtype=[('NuEnergy', float), 
                                 ('CosNuZenith', float), 
                                 ('PrimaryType', int)
                                ]
                         )
        mcf.close()
        return mcdata

    def _make_initial_data(self):
        pass

    def initialize_nuSQuIDS(self):
        czens, energies, initial_flux = self._make_initial_data()
        interactions = True
        self.nsq_atm = nsq.nuSQUIDSAtm(czens, energies*pc.GeV, 3, nsq.NeutrinoType.both, interactions)
        self.nsq_atm.Set_initial_state(initial_flux, nsq.Basis.flavor)
        
        self.nsq_atm.Set_rel_error(1.0e-15)
        self.nsq_atm.Set_abs_error(1.0e-15)
        self.nsq_atm.EvolveState()

    def get_flux(self, cz, e, ptype):
        pass
    
    def interp_mc(self):
        mcflux = [self.get_flux(*tup) for tup in zip(self.mcdata['CosNuZenith'],
                                                     self.mcdata['NuEnergy'],
                                                     self.mcdata['PrimaryType']
                                                    )
                 ]
        return mcflux

class SignalMCFluxMaker(BaseMCFluxMaker):

    def _make_initial_data(self):
        pg       = PathGen(self.mcpath)
        fluxfile = np.load('/data/user/jlazar/solar_WIMP/data/charon_fluxes/%s_1AU_BRW_dn_dz.npy' % (self.fluxtype))
            
        czens    = np.linspace(-1, 0.2, 100)
        energies = fluxfile['Energy']

        initial_flux = np.zeros((len(czens), len(energies), 2, 3), dtype=float)
        for ic in range(len(czens)):
            initial_flux[ic,:,0,0] = fluxfile['nu_e']
            initial_flux[ic,:,1,0] = fluxfile['nu_e_bar']
            initial_flux[ic,:,0,1] = fluxfile['nu_mu']
            initial_flux[ic,:,1,1] = fluxfile['nu_mu_bar']
            initial_flux[ic,:,0,2] = fluxfile['nu_tau']
            initial_flux[ic,:,1,2] = fluxfile['nu_tau_bar']
        del fluxfile
        return czens, energies, initial_flux

    def get_flux(self, cz, e, ptype):
        m = float(self.fluxtype.split('-')[1][1:])
        if e>m or e<10:
            return 0
        elif ptype==14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 0)
        elif ptype==-14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 1)
        else:
            print('wrong ptype doggo')
            return 0

class ConvNumuMCFluxMaker(BaseMCFluxMaker):

    def _make_initial_data(self):
        surface_flux_file = np.genfromtxt('/data/user/jlazar/solar_WIMP/data/AIRS_flux_sib_HG.dat')
        czens             = surface_flux_file[:,0][::350]
        energies          = surface_flux_file[:,1][:350]
        
        initial_flux = np.zeros((len(czens), len(energies), 2, 3), dtype=float)
        for ic in range(len(czens)):
            initial_flux[ic,:,0,0] = surface_flux_file[:,2][350*ic:350*(ic+1)]
            initial_flux[ic,:,1,0] = surface_flux_file[:,3][350*ic:350*(ic+1)]
            initial_flux[ic,:,0,1] = surface_flux_file[:,4][350*ic:350*(ic+1)]
            initial_flux[ic,:,1,1] = surface_flux_file[:,5][350*ic:350*(ic+1)]
            initial_flux[ic,:,0,2] = 0.0
            initial_flux[ic,:,1,2] = 0.0
        del surface_flux_file

        return czens, energies, initial_flux

    def get_flux(self, cz, e, ptype):
        if ptype==14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 0)
        elif ptype==-14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 1)
        else:
            print('wrong ptype doggo')
            return 0

class SolarAtmMCFluxMaker(BaseMCFluxMaker):
    
    def _make_initial_data(self):
        pp_HG_nu    = np.genfromtxt('/data/user/jlazar/solar_WIMP/data/solar_atm/PostPropagation/SIBYLL2.3_pp_HillasGaisser_H4a_nu.txt')
        pp_HG_nubar = np.genfromtxt('/data/user/jlazar/solar_WIMP/data/solar_atm/PostPropagation/SIBYLL2.3_pp_HillasGaisser_H4a_nubar.txt')
        czens       = np.linspace(-1, 0.2, 100)
        energies    = pp_HG_nu[:,0]
        initial_flux = np.zeros((len(czens), len(energies), 2, 3))

        for ic in range(len(czens)):
            initial_flux[ic,:,0,0] = pp_HG_nu[:,1]
            initial_flux[ic,:,1,0] = pp_HG_nu[:,2]
            initial_flux[ic,:,0,1] = pp_HG_nu[:,3]
            initial_flux[ic,:,1,1] = pp_HG_nubar[:,1]
            initial_flux[ic,:,0,2] = pp_HG_nubar[:,2]
            initial_flux[ic,:,1,2] = pp_HG_nubar[:,3]
        del pp_HG_nu
        del pp_HG_nubar
        return czens, energies, initial_flux

    def get_flux(self, cz, e, ptype):
        if e>1e5:
            return 0
        elif ptype==14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 0)
        elif ptype==-14:
            return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 1)
        else:
            print('wrong ptype doggo')
            return 0


if __name__=='__main__':
    args      = initialize_args()
    print('bong')
    if args.fluxtype=='conv-numu':
        fluxmaker = ConvNumuMCFluxMaker(args.mcfile, args.fluxtype)
    elif args.fluxtype=='solar-atm':
        fluxmaker = SolarAtmMCFluxMaker(args.mcfile, args.fluxtype)
    else:
        fluxmaker = SignalMCFluxMaker(args.mcfile, args.fluxtype)
    print('bang')
    fluxmaker.initialize_nuSQuIDS()
    print('bloop')
    mcflux    = fluxmaker.interp_mc()
    pg        = PathGen(args.mcfile)
    np.save(pg.get_mc_dn_dz_path(args.fluxtype)+'test', mcflux)
