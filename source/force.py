import numpy as np

from allvar import *


def distance(r1, r2):
    """Return Euclidean distance between positions"""
    return np.sqrt((np.array(r1) - np.array(r2))**2.)


def drdt(r, v):
    """Return position derivative

    :param r: shape: (x_earth, y_earth, x_jupiter, y_jupiter))
    :param v: shape: (vx_earth, vy_earth, vx_jupiter, vy_jupiter) 
    :return: velocities
    """
    return v

def dvdt(r, v, eps=1.e-32):
    """Return position derivative

    Central star have fixed position at (0, 0)

    :param r: shape: (x_earth, y_earth, x_jupiter, y_jupiter)
    :param v: shape: (vx_earth, vy_earth, vx_jupiter, vy_jupiter)
    :return: accelerations
    """
    # Normalization factors
    SOLAR_EARTH_NORM = G*SOLAR_MASS*EARTH_MASS*YEAR_TO_SEC**2./AU_TO_CM**3.
    SOLAR_JUPITER_NORM = G*SOLAR_MASS*JUPITER_MASS*YEAR_TO_SEC**2./AU_TO_CM**3.
    EARTH_JUPITER_NORM = G*EARTH_MASS*JUPITER_MASS*YEAR_TO_SEC**2./AU_TO_CM**3.

    # Geometric measurements
    r_se, r_sj, r_ej = r[:2], r[2:], r[2:] - r[:2]
    dist_se = distance((0, 0), r_se)
    dist_sj = distance((0, 0), r_sj)
    dist_ej = distance(r_se, r_sj)

    theta_se = np.math.athan(np.abs(r_se[1])/(np.abs(r_se[0]) + eps))
    theta_sj = np.math.athan(np.abs(r_sj[1])/(np.abs(r_sj[0]) + eps))
    theta_ej = np.math.athan(np.abs(r_ej[1])/(np.abs(r_ej[0]) + eps))
    
    # Unit force functionals 
    f_se = (
        -np.sign(r_se)
        * np.array(
            [
                np.cos(theta_se)/(dist_se + eps)**2.
                np.sin(theta_se)/(dist_se + eps)**2.
            ]
        )
    )
    f_sj = (
        -np.sign(r_sj)
        * np.array(
            [
                np.cos(theta_sj)/(dist_sj + eps)**2.
                np.sin(theta_sj)/(dist_sj + eps)**2.
            ]
        )
    )
    f_ej = (
        -np.sign(r_ej)
        * np.array(
            [
                np.cos(theta_ej)/(dist_ej + eps)**2.
                np.sin(theta_ej)/(dist_ej + eps)**2.
            ]
        )
    )

    # Accelearations
    a_e = SOLAR_EARTH_NORM*f_se - EARTH_JUPITER_NORM*f_ej 
    a_j = SOLAR_JUPITER_NORM*f_sj + EARTH_JUPITER_NORMf_ej

    return np.hstack([a_e, a_j])
