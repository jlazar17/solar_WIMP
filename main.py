from background_gamma import BackgroundGamma
from signal_gamma import SignalGamma

def main(mcpath, signal_flux, _skip=1):
    signal_gamma = SignalGamma(mcpath, signal_flux, _skip=_skip)
    signal_gamma.set_gamma_hist()
    signal_gamma.save_gamma_hist()
    bg_gamma = BackgroundGamma(mcpath, _skip=_skip)
    bg_gamma.set_gamma_hist()
    bg_gamma.save_gamma_hist()

if __name__=='__main__':
    from sys import argv as args
    main(args[1], args[2], _skip=1)
