# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
from typing import Callable

import numpy as np
import numpy.random as rd

from model.hydrogen import a_0

MONTE_CARLO_SAMPLING: int = 7_000


def uniform_spherical_interval(r_a: float = 0, r_b: float = a_0, theta_a: float = 0,
                               theta_b: float = np.pi * 2, phi_a: float = 0, phi_b: float = np.pi,
                               samples=MONTE_CARLO_SAMPLING):
    """
    :param r_a: initial radius
    :param r_b: final radius
    :param theta_a: initial polar angle
    :param theta_b: final polar angle
    :param phi_a: initial azimuthal angle
    :param phi_b: final azimuthal angle
    :param samples: number of elements given the ranges [r_a,r_b] [theta_a,theta_b] [phi_a,phi_b]
    :return: tree randomly uniform intervals
    """
    return rd.uniform(r_a, r_b, samples), rd.uniform(theta_a, theta_b, samples), rd.uniform(phi_a, phi_b, samples)


def spherical_interval(r_a: float = 0, r_b: float = a_0, theta_a: float = 0,
                       theta_b: float = np.pi * 2, phi_a: float = 0, phi_b: float = np.pi,
                       samples=MONTE_CARLO_SAMPLING):
    """
    :param r_a: initial radius
    :param r_b: final radius
    :param theta_a: initial polar angle
    :param theta_b: final polar angle
    :param phi_a: initial azimuthal angle
    :param phi_b: final azimuthal angle
    :param samples: number of elements given the ranges [r_a,r_b] [theta_a,theta_b] [phi_a,phi_b]
    :return: tree linearly uniform intervals
    """
    return np.linspace(r_a, r_b, samples), np.linspace(theta_a, theta_b, samples), np.linspace(phi_a, phi_b, samples)


def spherical_monte_carlo_integral(f: Callable[[float, float, float], float], r_a: float = 0, r_b: float = a_0,
                                   theta_a: float = 0, theta_b: float = np.pi * 2,
                                   phi_a: float = 0, phi_b: float = np.pi, samples=MONTE_CARLO_SAMPLING):
    """
    :param f: the function to be evaluated in spherical coordinates
    :param r_a: initial radius
    :param r_b: final radius
    :param theta_a: initial polar angle
    :param theta_b: final polar angle
    :param phi_a: initial azimuthal angle
    :param phi_b: final azimuthal angle
    :param samples: the number of samples in the integral
    :return: the approximate integral value for f
    """
    evaluated_volume: float = (((r_b - r_a) ** 3) / 3) * (theta_b - theta_a) * (-np.cos(phi_b) + np.cos(phi_a))
    r, theta, phi = uniform_spherical_interval(r_a, r_b, theta_a, theta_b, phi_a, phi_b, samples)
    summation: float = 0
    for r_i, theta_i, phi_i in zip(r, theta, phi):
        summation += f(r_i, theta_i, phi_i)
    return (evaluated_volume * summation) / samples, r, theta, phi, evaluated_volume


def spherical_monte_carlo_estimation_error(f: Callable[[float, float, float], float],
                                           r: np.array, theta: np.array, phi: np.array, volume: float,
                                           samples: int = MONTE_CARLO_SAMPLING):
    """
    :param f: The function evaluated in the spherical monte carlo integral
    :param r: the random values for the radius
    :param theta: the random values for the polar angle
    :param phi: the random values for the azimuthal angle
    :param volume: the spherical volume
    :param samples: the samples requested in the spherical monte carlo integral
    :return: The error estimation for the spherical monte carlo integral with previous parameters
    """
    summation: float = 0
    for r_i, theta_i, phi_i in zip(r, theta, phi):
        summation += f(r_i, theta_i, phi_i)
    summation /= samples

    sqrt_var_f: float = 0
    for r_i, theta_i, phi_i in zip(r, theta, phi):
        sqrt_var_f += (f(r_i, theta_i, phi_i) - summation) ** 2
    sqrt_var_f = np.sqrt(sqrt_var_f / samples - 1)
    
    return sqrt_var_f * volume / np.sqrt(samples)


def trapezoidal_rule(f: Callable[[float], float], f_second: Callable[[float], float],
                     a: float, b: float, n: int) -> float:
    """
    :param f: the function to evaluate
    :param f_second: the second derivative
    :param a: the lower integration limit
    :param b: the upper integration limit
    :param n: the number of sub-intervals
    :return: the approximate value for the integral of f
    """
    h = (b - a) / n
    xi = (b - a) / 2  # rd.uniform(a + 0.1, b - 0.1)
    summation = 0
    for x_i in np.linspace(start=a, stop=b - h, num=n):
        summation += f(x_i)
    return (h / 2) * (summation + (f(b) - f(a))) - (b - a) * ((h ** 2) / 12) * f_second(xi)


def midpoint_rule(f: Callable[[float], float], f_second: Callable[[float], float],
                  a: float, b: float, n: int) -> float:
    """
    :param f: the function to evaluate
    :param f_second: the second derivative
    :param a: the lower integration limit
    :param b: the upper integration limit
    :param n: the number of sub-intervals
    :return: the approximate value for the integral of f
    """
    h = (b - a) / n
    xi = (b - a) / 2  # rd.uniform(a + 0.1, b - 0.1)
    summation = 0
    for x_i in np.linspace(start=a, stop=b - h, num=n):
        summation += f(x_i + (h * 0.5))
    return h * summation + (b - a) * ((h ** 2) / 24) * f_second(xi)


def quadrature_rule(f: Callable[[float], float], f_second: Callable[[float], float],
                    a: float, b: float, n: int) -> float:
    """
    :param f: the function to evaluate
    :param f_second: the second derivative
    :param a: the lower integration limit
    :param b: the upper integration limit
    :param n: the number of sub-intervals
    :return: the approximate value for the integral of f
    """
    return ((2 * midpoint_rule(f, f_second, a, b, n)) + trapezoidal_rule(f, f_second, a, b, n)) / 3
