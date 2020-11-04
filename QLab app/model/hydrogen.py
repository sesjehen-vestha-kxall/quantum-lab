# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import numpy as np

a_0: float = 5.29177210903E-11  # First Bohr Radius


def hydrogen_1s(r: float, theta: float = 0, phi: float = 0):
    """
    :param r: the radius
    :param theta: the azimuthal angle
    :param phi: the polar angle
    :return: the value of psi_{1,0,0}(r, theta, phi)
    """
    return np.exp(-r / a_0) / (np.sqrt(np.pi) * a_0 ** (3 / 2))


def hydrogen_1s_second(r: float, theta: float = 0, phi: float = 0):
    """
    :param r: the radius
    :param theta: the azimuthal angle
    :param phi: the polar angle
    :return: the value of psi_{1,0,0}''(r, theta, phi)
    """
    return 4 * np.exp(-2*r / a_0) / (np.pi * (a_0**5))


def hydrogen_2s(r: float, theta: float = 0, phi: float = 0):
    """
    :param r: the radius
    :param theta: the azimuthal angle
    :param phi: the polar angle
    :return: the value of psi_{1,0,0}''(r, theta, phi)
    """
    return (np.exp(-r / a_0) * (2 - (r / a_0)) ** 2) / (32 * np.pi * (a_0 ** 3))


def hydrogen_2s_second(r: float, theta: float = 0, phi: float = 0):
    """
    :param r: the radius
    :param theta: the azimuthal angle
    :param phi: the polar angle
    :return: the value of psi_{1,0,0}''(r, theta, phi)
    """
    return (np.exp(-r / a_0) * (-8 * a_0 * r + 14 * a_0 ** 2 + r ** 2)) / (32 * np.pi * a_0 ** 7)
