import matplotlib as mpl
mpl.use("Agg")
import config
from physicsconstants import PhysicsConstants

param = PhysicsConstants()

print(config.NuFlux_Solar("Pythia",10, 10000, 200, "bb", 10000, param,location="Sunsfc"))
