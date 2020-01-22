import numpy as np
from numpy import cos as c
from numpy import sin as s


mc       = np.load("/data/user/jlazar/solar_WIMP/data/mcRecarray.npy")
nu_e     = mc["nuE"]
nu_zen   = mc["nuZen"]
nu_az    = mc["nuAz"]
reco_e   = mc["recoE"]
reco_zen = mc["recoZen"]
reco_az  = mc["recoAz"]

def meows_error(e):
    spl    = np.load("data/MEOWS_median_error_spline.npy").item()
    log10e = np.log10(e)
    return np.power(10, spl(log10e))

def opening_angle(zen1, az1, zen2, az2):
    return np.arccos(np.sin(zen1)*np.sin(zen2)*np.cos(az1-az2)+np.cos(zen1)*np.cos(zen2))

def rot1(nu_az, nu_zen, gen_az, gen_zen):
    new_x =  c(nu_zen,)*c(gen_az)*s(gen_zen) + s(nu_zen,)*c(gen_zen)
    new_y =  s(gen_az)*s(gen_zen)
    new_z = -s(nu_zen)*c(gen_az)*s(gen_zen) + c(nu_zen,)*c(gen_zen)
    new_az  = np.zeros(len(new_x))
    for i in range(len(new_x)):
        if new_x[i] > 0:
            new_az[i] = np.arctan(new_y[i] / new_x[i]) % (2*np.pi)
        else:
            new_az[i] = np.arctan(new_y[i] / new_x[i]) + np.pi
    new_zen = np.arccos(new_z)
    return new_az, new_zen

def rot2(nu_az, nu_zen, gen_az, gen_zen):
    new_x = c(nu_az)*c(gen_az)*s(gen_zen) - s(nu_az)*s(gen_az)*s(gen_zen)
    new_y = s(nu_az)*c(gen_az)*s(gen_zen) + c(nu_az)*s(gen_az)*s(gen_zen)
    new_z = c(gen_zen)
    new_az  = np.zeros(len(new_x))
    for i in range(len(new_x)):
        if new_x[i] > 0:
            new_az[i] = np.arctan(new_y[i] / new_x[i]) % (2*np.pi)
        else:
            new_az[i] = np.arctan(new_y[i] / new_x[i]) + np.pi
    new_zen = np.arccos(new_z)
    return new_az, new_zen

def rotate_coords(nu_az, nu_zen, gen_az, gen_zen):
    return rot2(nu_az, nu_zen, *rot1(nu_az,nu_zen, gen_az, gen_zen))

def gen_new_zen_az(scale):
    opening_angles  = opening_angle(nu_zen, nu_az, reco_zen, reco_az)
    gen_zen         = opening_angles * scale
    gen_az          = np.random.rand(len(opening_angles))*2*np.pi
    new_az, new_zen = rotate_coords(nu_az, nu_zen, gen_az, gen_zen)
    return new_az, new_zen
