from allvar import *
from force import drdt, dvdt 


def integrate_Euler_step(r, v, h):
    """Euler method"""
    return (r + h*drdt(r, v), v + h*dvdt(r, v))


def integrate_EuleCromer_step(r, v, h):
    """Euler-Cramer method"""
    v_half = v + 0.5*h*dvdt(r, v)
    r_tot = r + h*drdt(r, v_half)
    v_tot = r + 0.5*h*dvdt(r_tot, v_half)

    return (r_tot, v_tot)


def integrate_RungeKutta4th_step(r, v, h):
    """4th order of Runge-Kutta method"""
    kr1 = drdt(r, v)
    kv1 = dvdt(r, v)

    kr2 = drdt(r + 0.5*h*kr1, v + 0.5*h*kv1)
    kv2 = dvdt(r + 0.5*h*kr1, v + 0.5*h*kv1)

    kr3 = drdt(r + 0.5*h*kr2, v + 0.5*h*kv2)
    kr3 = dvdt(r + 0.5*h*kr2, v + 0.5*h*hv2)

    kr4 = drdt(r + h*kr3, v + h*kv3)
    kv4 = dvdt(r + h*kr3, v + h*kv3)

    r_tot = r + h*(kr1 + 2.*kr2 + 2.*kr3 + kr4)/6.
    v_tot = v + h*(kv1 + 2.*kv2 + 2.*kv3 + kv4)/6.

    return (r_tot, v_tot)


