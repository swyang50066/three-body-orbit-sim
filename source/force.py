import numpy as np

from allvar import *


def distance(r1, r2):
    """Return Euclidean distance between positions"""
    return np.sqrt(np.sum((np.array(r1) - np.array(r2))**2.))


def drdt(r, v):
    """Return position derivative

    :param r: shape: (x_earth, y_earth, x_jupiter, y_jupiter))
    :param v: shape: (vx_earth, vy_earth, vx_jupiter, vy_jupiter) 
    :return: velocities
    """
    return np.array(v)


def dvdt(r, v, eps=1.e-20):
    """Return position derivative

    Central star have fixed position at (0, 0)

    :param r: shape: (x_earth, y_earth, x_jupiter, y_jupiter)
    :param v: shape: (vx_earth, vy_earth, vx_jupiter, vy_jupiter)
    :return: accelerations
    """
    r, v = np.array(r), np.array(v)
    
    # Geometric measurements
    r_se, r_sj, r_ej = r[:2], r[2:], r[2:] - r[:2]
    dist_se = distance((0, 0), r_se)
    dist_sj = distance((0, 0), r_sj)
    dist_ej = distance(r_se, r_sj)

    theta_se = np.math.atan(np.abs(r_se[1])/(np.abs(r_se[0]) + eps))
    theta_sj = np.math.atan(np.abs(r_sj[1])/(np.abs(r_sj[0]) + eps))
    theta_ej = np.math.atan(np.abs(r_ej[1])/(np.abs(r_ej[0]) + eps))
    
    # Unit force functionals
    const_se = GG*(EARTH_MASS/SOLAR_MASS)
    f_se = -np.sign(r_se)*const_se*np.array(
        [      
            np.cos(theta_se)/(dist_se + eps)**2.,
            np.sin(theta_se)/(dist_se + eps)**2.
        ]
    )
    const_sj = GG*(JUPITER_MASS/SOLAR_MASS)
    f_sj = -np.sign(r_sj)*const_sj*np.array(
        [
            np.cos(theta_sj)/(dist_sj + eps)**2.,
            np.sin(theta_sj)/(dist_sj + eps)**2.
        ]
    )
    const_ej = GG*(EARTH_MASS*JUPITER_MASS/SOLAR_MASS**2.)
    f_ej = -np.sign(r_ej)*const_ej*np.array(
        [
            np.cos(theta_ej)/(dist_ej + eps)**2.,
            np.sin(theta_ej)/(dist_ej + eps)**2.
        ]
    )
    
    return np.hstack([
        (f_se - f_ej)/(EARTH_MASS/SOLAR_MASS),
        (f_sj + f_ej)/(JUPITER_MASS/SOLAR_MASS),
    ])
