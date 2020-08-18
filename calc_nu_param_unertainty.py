import numpy as np
import argparse
import emcee

import nuSQUIDSpy as nsq
import charon
from uncertainty_controls import ee, CC_central, CC_err_up, CC_err_down, qr_ch_dict 
from physicsconstants import PhysicsConstants
pc = PhysicsConstants() 

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
    parser.add_argument('--ordering',
                        type=str,
                        help='normal (no) or inverted (io) ordering'
                       )
    parser.add_argument('-n',
                        type=int,
                        help='number of things to simulate. Default is 100',
                        default=100
                       )
    parser.add_argument('--pythia',
                        type=bool,
                        default=False,
                        help='force using PYTHIA instead of BRW. BRW must be used below 100 GeV'
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

def make_emcee_params(ordering, ):
    if ordering=='no':
        mus    = np.asarray([33.82, 48.3, 8.61, 222, 7.39e-5, 2.523e-3])
        sigmas = np.array([(0.76, 0.78), (1.9,1.1), (0.13,0.13), (28,38), (0.20e-5,0.21e-5), (0.030e-3,0.032e-3)]).T
    elif ordering=='io':
        mus    = np.array([33.82, 48.6, 8.65, 285, 7.39e-5, -2.509e-3])
        sigmas = np.array([(0.76, 0.78), (1.5,1.1), (0.12,0.13), (26,24), (0.20e-5,0.21e-5), (0.030e-3,0.032e-3)]).T
    else:
        print('invalid ordering option. Must be "no" or "io"')
        quit()
    return mus, sigmas

def lnprob(theta, mus, sigma0s, sigma1s, ):
    vals = np.zeros(len(theta), dtype=float)
    tups = [(theta[i], mus[i], sigma0s[i], sigma1s[i]) for i in range(len(theta))]
    for i, tup in enumerate(tups):
        vals[i] = np.log(assymmetric_gaussian(*tup))
    return np.sum(vals[np.isfinite(vals)])

def make_osc_params(mus, sigmas, nwalkers=1000, nsteps=500):
    ndim = len(mus)
    pos = [mus + 1e-3*mus*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mus, sigmas[0], sigmas[1]))
    sampler.run_mcmc(pos, nsteps)
    chain  = sampler.chain
    osc_params = chain[:,-1,:]
    return osc_params

def calc_flux(ch, m, theta_12, theta_23, theta_13, delta_m_12, delta_m_13, delta, xsec, e_min, nodes, param, use_pythia):
    e_max = m
    dn_dz = np.zeros((2, nodes))
    if use_pythia:
        ws2qr = {5:'bb', 8:'WW', 11:'tautau'}
        flux_path = '/data/user/qliu/DM/DMFlux/Pythia/no_EW/secluded/Sun/results/%s_%d_Sun.dat' % (ws2qr[ch], m)
        flux = charon.NuFlux(qr_ch_dict[ch], m, nodes, Emin=e_min, Emax=e_max, process='ann',
                             theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                             delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                             interactions=True, xsec=xsec, bins=200, pathFlux=flux_path
                            )
    else:
        flux = charon.NuFlux(qr_ch_dict[ch], m, nodes, Emin=e_min, Emax=e_max, process='ann',
                             theta_12=theta_12, theta_13=theta_13, theta_23=theta_23,
                             delta=delta, delta_m_12=delta_m_12, delta_m_13=delta_m_13,
                             interactions=True, xsec=xsec, bins=200
                            )
    sun_flux   = flux.Sun('SunSurface', zenith=0., avg=True)
                           
    return sun_flux

def main(ch, m, ordering, nwalkers, param, use_pythia):

    savedir = '/data/user/jlazar/solar_WIMP/data/param_uncertainties/'    

    mus, sigmas = make_emcee_params(ordering)
    
    osc_params = make_osc_params(mus, sigmas, nwalkers=nwalkers)
    
    th12, th23, th13, delta, m12, m13 = tuple(mus)
    xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP_py2/sources/nuSQuIDS/data/xsections/csms.h5', ee, [1.])
    nominal_sun = calc_flux(ch, m, th12, th23, th13, m12, m13, delta, xsec, 10, 200, param, use_pythia)
    
    if use_pythia:
        np.save("%s/%d_%d_%s_nominal_pythia.npy" % (savedir, ch, m, ordering), nominal_sun)
    else:
        print('blonk')
        np.save("%s/%d_%d_%s_nominal.npy" % (savedir, ch, m, ordering), nominal_sun)


    for osc_param in osc_params:
        th12, th23, th13, delta, m12, m13 = tuple(osc_param)
        xsec=nsq.ScaledNeutrinoCrossSections('/home/jlazar/programs/GOLEM_SOLAR_WIMP_py2/sources/nuSQuIDS/data/xsections/csms.h5', ee, [1.])
        print(use_pythia)
        if use_pythia:
            savefname = '%d_%d_%s_%f_%f_%f_%f_%f_%f_pythia.npy' % (ch, m, ordering, th12, th23, th13, m12, m13, delta)
        else:
            savefname = '%d_%d_%s_%f_%f_%f_%f_%f_%f.npy' % (ch, m, ordering, th12, th23, th13, m12, m13, delta)
            print(savefname)
        sun_flux = calc_flux(ch, m, th12, th23, th13, m12, m13, delta, xsec, 10, 200, param,use_pythia)
        np.save("%s/%s" % (savedir, savefname), sun_flux)


if __name__=='__main__':
    args = initialize_parser()
    ch       = args.ch
    m        = args.m
    ordering = args.ordering
    n        = args.n
    use_pythia = args.pythia
    print(use_pythia)
    seed = hash(str(ch)+str(m)+ordering) % 2**32
    np.random.seed(seed)
    main(ch, m, ordering, n, pc, use_pythia)
