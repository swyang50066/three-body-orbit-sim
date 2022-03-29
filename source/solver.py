from force import drdt, dvdt 


def integrate_Euler_step(r, v, h):
    """Euler method"""
    return (r + h*drdt(r, v), v + h*dvdt(r, v))


def integrate_EulerCromer_step(r, v, h):
    """Euler-Cramer method"""
    # Integrate 
    v_half = v + 0.5*h*dvdt(r, v)
    r_tot = r + h*drdt(r, v_half)
    v_tot = r + 0.5*h*dvdt(r_tot, v_half)

    return (r_tot, v_tot)


def integrate_RungeKutta4th_step(r, v, h):
    """4th order of Runge-Kutta method"""
    # Runge-Kutta coefficients
    r_k1 = drdt(r, v)
    v_k1 = dvdt(r, v)

    r_k2 = drdt(r + 0.5*h*r_k1, v + 0.5*h*v_k1)
    v_k2 = dvdt(r + 0.5*h*r_k1, v + 0.5*h*v_k1)

    r_k3 = drdt(r + 0.5*h*r_k2, v + 0.5*h*v_k2)
    v_k3 = dvdt(r + 0.5*h*r_k2, v + 0.5*h*v_k2)

    r_k4 = drdt(r + h*r_k3, v + h*v_k3)
    v_k4 = dvdt(r + h*r_k3, v + h*v_k3)

    # Updates
    r_tot = r + h*(r_k1 + 2.*r_k2 + 2.*r_k3 + r_k4)/6.
    v_tot = v + h*(v_k1 + 2.*v_k2 + 2.*v_k3 + v_k4)/6.

    return (r_tot, v_tot)


