import numpy as np

def logfu_to_ppb(logfu):
    """
    Convert log10(fraction unbound) back to % PPB.
    """

    fu = 10 ** logfu

    ppb = (1.0 - fu) * 100.0

    return np.clip(ppb, 0.0, 100.0)


def log_to_vd(log_vd):
    """
    Inverse log transform for Vd.
    """

    return np.exp(log_vd)


def logit_to_bioavailability(logit):
    """
    Convert logit(F) back to bioavailability %.
    """

    frac = 1.0 / (1.0 + np.exp(-logit))

    return frac * 100.0