import file_gen as fg

thing = fg.File_Gen("/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/HoleIce1/Ares/IC86.AVG/Merged/Ares_IC86.AVG_1.27_lite_platinum_99510.h5")

print(thing.mcfile)
print(thing.mctype)
print(thing.get_mcname())
print(thing.get_mc_dn_dz_path("ch8-m1000")) 
