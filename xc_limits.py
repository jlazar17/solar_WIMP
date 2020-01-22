scDict = {5 :(3.06e-03, 2.82e-06, 2.59e-03, 2.00e-06),
          8 :(3.76e-05, 3.49e-08, 6.80e-05, 5.28e-08),
          11:(1.46e-05, 1.35e-08, 2.07e-05, 1.60e-08)
         }

class Limits():

    def __init__(self, ch):
        
        self.SD500,self.SI500,self.SD1000,self.SI1000 = scDict[ch]
