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

def pic50_to_ic50(pic50):
    """
    Convert pIC50 to IC50 in nM.
    """
    pic50 = np.asarray(
        pic50,
        dtype=np.float64,
    )

    return np.power(
        10.0,
        9.0 - pic50,
    )


def logfu_to_ppb(logfu):
    """
    Convert log10 fraction unbound to percentage
    plasma protein bound.
    """
    logfu = np.asarray(
        logfu,
        dtype=np.float64,
    )

    fu = np.power(
        10.0,
        logfu,
    )

    return 100.0 * (
        1.0 - fu
    )


def log_to_vd(log_vd):
    """
    Reverse the natural-log Vd transform.
    """
    log_vd = np.asarray(
        log_vd,
        dtype=np.float64,
    )

    return np.exp(
        log_vd
    )


def logit_to_bioavailability(logit_f):
    """
    Convert logit-scaled bioavailability back to
    bioavailability percentage.
    """
    logit_f = np.asarray(
        logit_f,
        dtype=np.float64,
    )

    fraction = np.empty_like(
        logit_f,
        dtype=np.float64,
    )

    positive = logit_f >= 0

    fraction[positive] = (
        1.0 /
        (
            1.0 +
            np.exp(
                -logit_f[positive]
            )
        )
    )

    exp_values = np.exp(
        logit_f[~positive]
    )

    fraction[~positive] = (
        exp_values /
        (
            1.0 +
            exp_values
        )
    )

    return 100.0 * fraction