import numpy as np

def openingAngle(zen1, az1, zen2, az2):
    return np.arccos(np.sin(zen1)*np.sin(zen2)*np.cos(az1-az2)+np.cos(zen1)*np.cos(zen2))

class MonteCarlo():

    def __init__(self, arrPath, nuType):
        assert (nuType in ["nu", "nuBar"]), "Invalid nuType"
        print("papp")
        self.recArray = np.load(arrPath)
        print("peep")
        self.nuType   = nuType
        print("piip")
        self.setIndices()
        self.loadH5Info()

    def setIndices(self):
        if self.nuType=="nu":
            self.i = np.where(self.recArray["i"]==14)[0]
        elif self.nuType=="nuBar":
            self.i = np.where(self.recArray["i"]==-14)[0]

    def loadH5Info(self):
        self.nuZen     = self.recArray["nuZen"][self.i]
        self.nuAz      = self.recArray["nuAz"][self.i]
        self.nuE       = self.recArray["nuE"][self.i]
        self.recoZen   = self.recArray["recoZen"][self.i]
        self.recoAz    = self.recArray["recoAz"][self.i]
        self.recoE     = self.recArray["recoE"][self.i]
        self.oneWeight = self.recArray["oneWeight"][self.i] * 1.e-4

    def setTrueGamma(self, zen, az):
        self.trueGamma = openingAngle(zen, az, self.nuZen, self.nuAz)

    def setRecoGamma(self, zen, az):
        self.recoGamma = openingAngle(zen, az, self.recoZen, self.recoAz)

    def setNuZen(self, nuZen):
        self.nuZen = nuZen

    def setNuAz(self, nuAz):
        self.nuAz = nuAz

    def setNuE(self, nuE):
        self.nuE = nuE

    def setRecoZen(self, recoZen):
        self.recoZen = recoZen

    def setRecoAz(self, recoAz):
        self.recoAz = recoAz

    def setRecoE(self, recoE):
        self.recoE = recoE

    def setOneWeight(self, oneWeight):
        self.oneWeight = oneWeight
