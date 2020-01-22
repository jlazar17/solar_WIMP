import argparse

parser = parser = argparse.ArgumentParser()
parser.add_argument("-f",
                    type=float,
                    default=1.,
                    help="factor to rescale MEOWS"
                   )
parser.add_argument("--ch",
                    type=int,
                    help="WIMPSim channel number bb:5, WW:8, tautau:11"
                   )
parser.add_argument("-m",
                    type=int,
                    help="Dark matter mass"
                   )
parser.add_argument("-n",
                    type=int,
                    help="run number"
                   )
parser.add_argument("--nt",
                    type=str,
                    help="neutrino type (nu or nuBar)"
                   )
args   = parser.parse_args()
factor = args.f
ch     = args.ch
m      = args.m
nt     = args.nt
n      = args.n

print(factor)
print(ch)
print(m)
print(nt)
print(n)
