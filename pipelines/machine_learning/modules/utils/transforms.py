import numpy as np

def ic50_to_pic50(x):
    return -np.log10(x * 1e-9)  # assume values are in nM

def ppb_to_logfu(ppb_percent):
    """
    Convert % plasma protein bound to log10(fraction unbound).

    Example:
        99% bound -> fu=0.01 -> log10(fu)=-2
    """

    fu = 1.0 - (ppb_percent / 100.0)

    # numerical stability
    fu = np.clip(fu, 1e-6, 1.0)

    return np.log10(fu)


def vd_to_log(vd):
    """
    Log-transform volume of distribution.
    Assumes positive values.
    """

    vd = np.clip(vd, 1e-6, None)

    return np.log(vd)


def bioavailability_to_logit(f_percent):
    """
    Convert bioavailability % to logit space.
    """

    frac = f_percent / 100.0

    # avoid infs
    frac = np.clip(frac, 1e-4, 1 - 1e-4)

    return np.log(frac / (1.0 - frac))