import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import os
import argparse
import nuSQUIDSpy as nsq
import emcee
from scipy.interpolate import splrep, splev, interp1d

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

def assymmetric_gaussian(x, mu, sigma0, sigma1):
    norm = 2 / (np.sqrt(np.pi)*(sigma0+sigma1))
    if x <= mu:
        val = norm * np.exp(-np.power(x-mu,2)/np.power(sigma0, 2))
    else:
        val = norm * np.exp(-np.power(x-mu,2)/np.power(sigma1, 2))
    return val

def make_emcee_params():
    mus    = np.ones(6)
    sigmas = np.array([(CC_err_up[i]-1,1-CC_err_down[i]) for i in range(6)]).T
    print(sigmas)
    return mus, sigmas

def lnprob(theta, mus, sigma0s, sigma1s, ):
    vals = np.zeros(len(theta), dtype=float)
    tups = [(theta[i], mus[i], sigma0s[i], sigma1s[i]) for i in range(len(theta))]
    for i, tup in enumerate(tups):
        vals[i] = np.log(assymmetric_gaussian(*tup))
    return np.sum(vals[np.isfinite(vals)])

def make_scale_params(mus, sigmas, nwalkers=1000, nsteps=500):
    ndim = len(mus)
    pos = [mus + 1e-3*mus*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mus, sigmas[0], sigmas[1]))
    sampler.run_mcmc(pos, nsteps)
    chain  = sampler.chain
    osc_params = chain[:,-1,:]
    return np.where(osc_params>=0, osc_params, 0)

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
    return dn_dz

if __name__=='__main__':
    args       = initialize_parser()
    e_min      = 10 # GeV
    ch         = args.ch
    m          = args.m 
    xs_model  = args.xs
    ordering   = args.ordering
    theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta = set_osc_params(ordering)
    mus, sigmas = make_emcee_params()
    nodes      = 200
    seed = (hash(str(ch)+str(m)+ordering)+1) % 2**32
    np.random.seed(seed)
    if xs_model=='error':
        ee = ee[:6]
        xs_scales = make_scale_params(mus, sigmas)
        for xs_scale in xs_scales:
            ee_fine = np.linspace(ee[0], ee[-1], 100)
            tck = splrep(ee, xs_scale, k=1)
            f = interp1d(ee, xs_scale)
            xs_scale_fine = splev(ee_fine, tck)
            print(xs_scale_fine)
            xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee_fine*1e9, xs_scale_fine)
            dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
            np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]), dn_dz)
        xs_scales = make_scale_arrs(CC_err_up)
       # for i, xs_scale in enumerate(xs_scales):
       #     xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee*1e9, xs_scale)
       #     dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
       #     print('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]))
       #     np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scale[0]), dn_dz)
        quit()
    elif xs_model=='csms':
        xsec = '/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5'
    elif xs_model=='nusigma':
        xsec='/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/nusigma_'
    elif xs_model=='nominal':
        xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP/sources/nuSQuIDS/data/xsections/csms.h5', ee, [1])
        dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
        np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_nominal.npy' % (ch,m,ordering,), dn_dz)
        quit()
    else:
        print('invalid xs model. Only "nusigma" and "csms" supported at this time')
        quit()
    dn_dz = calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, e_min, nodes, param, xsec)
    np.save('/data/user/jlazar/solar_WIMP/data/xs_uncertainties/%d_%d_%s_%f.npy' % (ch,m,ordering,xs_scales[0]), dn_dz)
