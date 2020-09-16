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

class MCFluxMaker:

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
        if self.fluxtype=='conv-numu':
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
        else:
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
        return czens, energies, initial_flux

    def initialize_nuSQuIDS(self):
        czens, energies, initial_flux = self._make_initial_data()
        interactions = True
        self.nsq_atm = nsq.nuSQUIDSAtm(czens, energies*pc.GeV, 3, nsq.NeutrinoType.both, interactions)
        self.nsq_atm.Set_initial_state(initial_flux, nsq.Basis.flavor)
        
        self.nsq_atm.Set_rel_error(1.0e-15)
        self.nsq_atm.Set_abs_error(1.0e-15)
        self.nsq_atm.EvolveState()

    def get_flux(self, cz, e, ptype):
        if self.fluxtype=='conv-numu':
            if ptype==14:
                return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 0)
            elif ptype==-14:
                return self.nsq_atm.EvalFlavor(1, cz, e*pc.GeV, 1)
            else:
                print('wrong ptype doggo')
                return 0
        else:
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
    
    def interp_mc(self):
        mcflux = [self.get_flux(*tup) for tup in zip(self.mcdata['CosNuZenith'],
                                                     self.mcdata['NuEnergy'],
                                                     self.mcdata['PrimaryType']
                                                    )
                 ]
        return mcflux

if __name__=='__main__':
    args      = initialize_args()
    print('bong')
    fluxmaker = MCFluxMaker(args.mcfile, args.fluxtype)
    print('bang')
    fluxmaker.initialize_nuSQuIDS()
    print('bloop')
    mcflux    = fluxmaker.interp_mc()
    pg        = PathGen(args.mcfile)
    np.save(pg.get_mc_dn_dz_path(args.fluxtype), mcflux)
