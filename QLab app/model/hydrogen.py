# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import numpy as np
import numpy.random as rd

MONTECARLO_SAMPLING = 60
N = MONTECARLO_SAMPLING ** 3

a_0 = 5.29177210903e-11 #First Bohr Radious

def hydrogen_1s(r: float, theta: float=0, phi: float=0):
    return np.exp(-r/a_0) / (np.sqrt(np.pi) * a_0**(3/2))

def random_interval_values(r_a: float=0, r_b: float=a_0, theta_a: float=0,
        theta_b: float=np.pi*2, phi_a: float=0, phi_b: float=np.pi*2, samples=MONTECARLO_SAMPLING):
    # Selected intervals
    return rd.uniform(r_a, r_b, samples),\
            rd.uniform(theta_a, theta_b, samples),\
            rd.uniform(phi_a, phi_b, samples)

def linear_interval_values(r_a: float=0, r_b: float=a_0, theta_a: float=0,
        theta_b: float=np.pi*2, phi_a: float=0, phi_b: float=np.pi*2, samples=MONTECARLO_SAMPLING):
    # Selected intervals
    return np.linspace(r_a, r_b, samples),\
            np.linspace(theta_a, theta_b, samples),\
            np.linspace(phi_a, phi_b, samples)

def probability(psi, r_a: float=0, r_b: float=a_0, theta_a: float=0,
        theta_b: float=np.pi*2, phi_a: float=0, phi_b: float=np.pi*2):
    evaluated_volume: float = 4 * np.pi * ((r_b - r_a) ** 3) / 3
    r, theta, phi = random_interval_values(r_a, r_b, theta_a, theta_b, phi_a, phi_b)
    summation: float = 0
    for r_i in r:
        for theta_i in theta:
            for phi_i in phi:
                summation += np.abs(psi(r_i, theta_i, phi_i)) ** 2
    return (evaluated_volume * summation) / N