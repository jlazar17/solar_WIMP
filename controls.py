import numpy as np
r_sun        = 6.9e10  # radius of sun in cm
au2cm = 1.496e13
delta_t = 60. * 30. # sec
jd_start = 2455349.5
jd_end = 2457906.5
n_days = jd_end - jd_start
n_step = delta_t / (3600. * 24.)
jds = np.linspace(jd_start, jd_end, int(n_days / n_step) + 1)
azimuths = np.random.rand(len(jds)) * 2 * np.pi # solar azimuth positions
gamma_bins = np.linspace(0, np.pi, 361) # bins of width one half degree
bins_per_decade = 10
high_exp       = 6.5
low_exp        = 0.5
e_bins         = np.logspace(low_exp, high_exp, int(bins_per_decade * (high_exp-low_exp) + 1))

chs = [5,8,11]
ms = [300,500,800,1000,2000,3000,10000]
mcfiles = [
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce3/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_97840.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce4/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_95655.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_95225.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce3/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14995.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce4/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99810.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce2/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_90260.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_37525.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce2/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Ares/IC86.AVG/Merged/Ares_IC86.AVG_0.97_lite_platinum_98000.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_Noise/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14615.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/Noise/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_90000.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_spice_mie/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.30_lite_platinum_18090.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/spice_mie/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.30_lite_platinum_73155.h5',
]

syst_dict = {
            'HoleIce1':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
                       ],
            'HoleIce3':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce3/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_97840.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce3/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14995.h5',
                       ],
            'HoleIce4':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce4/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_95655.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce4/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
                       ],
            'HoleIce5':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_95225.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_HoleIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
                       ],
            'BulkIce1':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99810.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
                       ],
            'BulkIce2':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce2/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_90260.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce2/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14990.h5',
                       ],
            'BulkIce5':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/BulkIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_37525.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_BulkIce5/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_15000.h5',
                       ],
            'Nominal':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Ares/IC86.AVG/Merged/Ares_IC86.AVG_0.97_lite_platinum_98000.h5',
                      ],
            'Noise':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_Noise/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_14615.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/Noise/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_90000.h5',
                    ],
            'spice_mie':[
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/LE_spice_mie/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.30_lite_platinum_18090.h5',
'/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/spice_mie/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.30_lite_platinum_73155.h5',
                        ],
           }

si_ann_rate_dict = {200  : 1.851394541398883e+28,
                    300  : 1.0262342643319446e+28,
                    400  : 6.546934291020078e+27,
                    500  : 4.5418315143731186e+27,
                    600  : 3.334488508098688e+27,
                    700  : 2.551030591392947e+27,
                    800  : 2.0139313842344916e+27,
                    900  : 1.6298244822123644e+27,
                    1000 : 1.3457320571730637e+27,
                    2000 : 3.6648784426051484e+26,
                    3000 : 1.674577361367118e+26,
                    4000 : 9.548136632702172e+25,
                    6000 : 4.3005454412986445e+25,
                    7000 : 3.171480260094724e+25,
                    8000 : 2.434979863489246e+25,
                    9000 : 1.9281149572541875e+25,
                    10000: 1.5644775205377196e+25
                   } # 1/sec (xs=1 pb)

sd_ann_rate_dict = {200  : 4.114212656452698e+25,
                    300  : 1.8483123974320923e+25,
                    400  : 1.044156901744328e+25,
                    500  : 6.697084065971152e+24,
                    600  : 4.6565477201402684e+24,
                    700  : 3.423814217222649e+24,
                    800  : 2.622728847911479e+24,
                    900  : 2.0730384857875426e+24,
                    1000 : 1.679607455216801e+24,
                    2000 : 4.202600488405267e+23,
                    3000 : 1.8680716988310504e+23,
                    4000 : 1.050822229370311e+23,
                    5000 : 6.7252930524968645e+22,
                    6000 : 4.6703267898238e+22,
                    7000 : 3.431240497707627e+22,
                    8000 : 2.627026669461829e+22,
                    9000 : 2.0756623642861236e+22,
                    10000: 1.681276619465538e+22
                   } # 1/sec (xs=1 pb)
