import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import os
import argparse
import nuSQUIDSpy as nsq

import charon
from physicsconstants import PhysicsConstants
from uncertainty_controls import ee, CC_central, CC_err_up, CC_err_down

param = PhysicsConstants()
qr_ch_dict = {5:"bb", 8:"WW", 11:"tautau"}
ws_cn_dict = {"bb":5, "WW":8, "tautau":11}
xsec_path = '/home/jlazar/program/GOLEMSOLARWIMP/sources/nuSQuIDS/resource/xsections/'


def initialize_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ch",
                        type=int,
                        help="WIMPSim channel number bb:5, WW:8, tautau:11"
                       )
    parser.add_argument("-m",
                        type=int,
                        help="Dark matter mass"
                       )
    parser.add_argument('--xs',
                        type=str,
                        help='cross section type ("nusigma" or "csms")'
                       )
    parser.add_argument('--ordering',
                        type=str,
                        help='normal (no) or inverted (io) ordering'
                       )
    args = parser.parse_args()
    return args

def make_scale_arrs(scales, n=100):
    xs_scale_factors = []
    std_devs = np.abs(1-np.asarray(scales))
    for std_dev in std_devs:
        xs_scale_factors.append(np.abs(np.random.normal(scale=std_dev, size=n)))
    xs_scale_factors = np.array(xs_scale_factors).T
    if scales[0]>1:
        return 1+xs_scale_factors
    else:
        return 1-xs_scale_factors
    
        

def set_osc_params(ordering):
    if ordering=='no':
        theta_12   = 33.82
        theta_23   = 48.3
        theta_13   = 8.61
        delta_m_12 = 0.0000739
        delta_m_13 = 0.002523
        delta      = 222.
    elif ordering=='io':
        theta_12   = 33.8200000
        theta_23   = 48.6
        theta_13   = 8.65
        delta_m_12 = 0.0000739
        delta_m_13 = -0.002509
        delta      = 285.0
    else:
        print('ordering %s not valid. Must be "no" or "io"' % ordering)
        quit()
    return (theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta)

def calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, params, xsec):
    e_max = m
    dn_dz = np.zeros((2, nodes))
    flux = charon.NuFlux(qr_ch_dict[ch], m, nodes, Emin=e_min, Emax=e_max, process='ann',
                         theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                         delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                         interactions=True, xsec=xsec, bins=200
                        )
    sun_flux   = flux.Sun('SunSurface', zenith=0.)
    dn_dz[0][:] = sun_flux['nu_mu']
    dn_dz[1][:] = sun_flux['nu_mu_bar']
#    f = config.NuFlux("Pythia", e_min, e_max, nodes, qr_ch_dict[ch], m, param,
#                           theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
#                           delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
#                           interactions=True, xsec=xsec_path, location='Sunsfc')
#    nu_mu_dn_dz     = np.asarray([tup[2] for tup in f]) * float(m)
#    nu_mu_bar_dn_dz = np.asarray([tup[5] for tup in f]) * float(m)
#    dn_dz[0][:] = nu_mu_dn_dz
#    dn_dz[1][:] = nu_mu_bar_dn_dz
    return dn_dz

if __name__=='__main__':
    args       = initialize_parser()
    e_min      = 10 # GeV
    ch         = args.ch
    m          = args.m 
    xs_model  = args.xs
    ordering   = args.ordering
    theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta = set_osc_params(ordering)
    nodes      = 200
    seed = (hash(str(ch)+str(m)+ordering)+1) % 2**32
    np.random.seed(seed)
    if xs_model=='error':
        xs_scales = make_scale_arrs(CC_err_down)
        for i, xs_scale in enumerate(xs_scales):
            print(xs_scale)
            xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee*1e9, xs_scale)
            dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
            print('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]))
            np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]), dn_dz)
        xs_scales = make_scale_arrs(CC_err_up)
       # for i, xs_scale in enumerate(xs_scales):
       #     xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee*1e9, xs_scale)
       #     dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
       #     print('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]))
       #     np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]), dn_dz)
        quit()
    elif xs_model=='csms':
        print('csms')
        xsec = '/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5'
    elif xs_model=='nusigma':
        print('nusigma')
        xsec='/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/nusigma_'
    elif xs_model=='nominal':
        xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee, [1])
        print('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_nominal.npy') % (ch,m,ordering,)
        dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
        np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_nominal.npy' % (ch,m,ordering,), dn_dz)
        quit()
    else:
        print('invalid xs model. Only "nusigma" and "csms" supported at this time')
        quit()
    dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
    np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scales[0]), dn_dz)
