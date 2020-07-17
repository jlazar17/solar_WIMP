from calc_sens import SensCalc
import numpy as np

#mc_files = [
#'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5',
#'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5']
mc_files = ['/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Ares/IC86.AVG/Merged/Ares_IC86.AVG_0.97_lite_platinum_98000.h5']
sc = SensCalc('ch8-m1000')
for f in mc_files:
    print('weep')
    sc.add_mub(f)
    print('woop')
    sc.add_mus(f)
sc.set_sens()
print(sc.sens)
