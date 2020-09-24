# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import numpy as np
import numpy.random as rd

MONTECARLO_VARIABLE_SAMPLING: int = 90

a_0 = 0.0529E-11 #First Bohr Radious

def maps_dict_append(d: dict, k, v):
    if k not in d:
        d[k] = {}
    d[k].append(v)

def hydrogen_1s(r: float, theta: float=0, phi: float=0):
    return np.exp(-r/a_0) / (np.sqrt(np.pi) * a_0**(3/2))

def random_coordinates(r_a: float=0, r_b: float=a_0, theta_a: float=0, theta_b: float=np.pi*2,
        phi_a: float=0, phi_b: float=np.pi*2, samples=MONTECARLO_VARIABLE_SAMPLING) -> set:
    # Selected intervals
    coords: set = { (0,0,0) }
    coords_number: int = samples**3
    while len(coords) < coords_number:
        coords.add((
            rd.uniform(r_a, r_b),
            rd.uniform(theta_a, theta_b),
            rd.uniform(phi_a, phi_b)
        ))
    return coords

def linear_interval_values(r_a: float=0, r_b: float=a_0, theta_a: float=0,
        theta_b: float=np.pi*2, phi_a: float=0, phi_b: float=np.pi*2, v_samples=MONTECARLO_VARIABLE_SAMPLING):
    # Selected intervals
    return np.linspace(r_a, r_b, v_samples),\
            np.linspace(theta_a, theta_b, v_samples),\
            np.linspace(phi_a, phi_b, v_samples)

def probability(psi, r_a: float=0, r_b: float=a_0, theta_a: float=0,
        theta_b: float=np.pi*2, phi_a: float=0, phi_b: float=np.pi*2):
    evaluated_volume: float = 4 * np.pi * ((r_b - r_a) ** 3)
    f_coords = random_coordinates(r_a, r_b, theta_a, theta_b, phi_a, phi_b)
    summation: float = 0
    for x_bar in f_coords:
        r_i, theta_i, phi_i = x_bar
        summation += np.abs(psi(r_i, theta_i, phi_i)) ** 2
    return (evaluated_volume * summation) / (MONTECARLO_VARIABLE_SAMPLING ** 3)