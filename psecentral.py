"""
psecentral.py

Utility functions for determining whether a given gamma value lies
within specific geometric or physical bounds.

"""
import pconstants

def central(gamma: float) -> bool:
    """
    Return True if abs(gamma) is within the central limit.

    Parameters
    ----------
    gamma : float
        The gamma value to test.

    Returns
    -------
    bool
        True if abs(gamma) <= CENTRAL_LIMIT, otherwise False.
    """
    return abs(gamma) <= pconstants.CENTRAL_GAMMA_LIMIT


def umbra_touch(gamma: float) -> bool:
    """
    Return True if abs(gamma) indicates umbra contact.

    Notes
    -----
    Currently this uses the same limit as `central()`. If these
    functions are intended to represent different physical conditions,
    adjust the limit accordingly.

    Parameters
    ----------
    gamma : float
        The gamma value to test.

    Returns
    -------
    bool
        True if abs(gamma) <= CENTRAL_LIMIT, otherwise False.
    """
    return abs(gamma) <= pconstants.CENTRAL_GAMMA_LIMIT


def exist(gamma: float) -> bool:
    """
    Return True if abs(gamma) is within the existence limit.

    Parameters
    ----------
    gamma : float
        The gamma value to test.

    Returns
    -------
    bool
        True if abs(gamma) <= EXIST_LIMIT, otherwise False.
    """
    return abs(gamma) <= pconstants.EXIST_GAMMA_LIMIT
