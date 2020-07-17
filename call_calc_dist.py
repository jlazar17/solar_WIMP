from recalc_energy_delta_theta_hist_bg import BackgroundGamma, SignalGamma
import numpy as np

sigG = SignalGamma('/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5', 'ch8-m1000', _skip=1000)
hist = sigG.gamma_calc()
np.save("/data/user/jlazar/solar_WIMP/recalc_sig_test.npy", hist)
bgG = BackgroundGamma('/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5', _skip=1000)
hist = bgG.gamma_calc()
np.save("/data/user/jlazar/solar_WIMP/recalc_bg_test.npy", hist)
