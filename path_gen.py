class PathGen():

    def __init__(self, mcpath):
        
        self.mcpath = mcpath
        self.mctype = mcpath.split("/")[7]

    def get_mcname(self):
        if self.mctype=="Systematics":
            mcname = self.mcpath.split("/")[8]
        elif self.mctype=="Nominal":
            mcname = "Nominal"
        else:
            print("mctype==%s is not a recognized mctype" % self.mctype)
            mcname=None
        return mcname

    def get_licfile_path(self):
        if self.mctype=="Nominal":
            licfile = "/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Gen/00001-01000/Generation_data.lic"
        elif self.mctype=="Systematics":
            licfile = "/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/%s/Gen/00001-01000/Generation_data.lic" % self.get_mcname()
        else:
            print("mctype==%s is not a recognized mctype" % self.mctype)
            licfile=None
        return licfile

    def get_mc_dn_dz_path(self, flux):
        return "/data/user/jlazar/solar_WIMP/data/mc_dn_dz/%s_%s_dn_dz.npy" % (flux, self.get_mcname())

    def get_ow_path(self):
        path = '/data/user/jlazar/solar_WIMP/data/mc_oneweight/%s.npy' % self.get_mcname()
        return path
 
    def get_mcpath(self):
        return self.mcpath

    def get_e_d_theta_path(self, fluxtype):
        return '/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/%s_%s_e_d_theta' % (fluxtype, self.get_mcname())

    def get_uninterp_flux_path(self, fluxtype):
        if fluxtype=='conv-numu':
            path = "/data/user/jlazar/solar_WIMP/data/AIRS_flux_sib_HG.dat"
        else:
            path =  '/data/user/jlazar/solar_WIMP/data/qr_dn_dz/%s_dn_dz.npy' % (fluxtype)
        return path
