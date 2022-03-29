import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from allvar import *
from solver import (
    integrate_Euler_step,
    integrate_EulerCromer_step,
    integrate_RungeKutta4th_step,
)


class Planet(object):
    """class object of planet state"""
    def __init__(self, r, v):
        self._pos = list([r])
        self._vel = list([v])

    @property
    def pos(self):
        return np.array(self._pos)

    @property
    def vel(self):
        return np.array(self._vel)
   
    @pos.setter
    def pos(self, Ipos):
        self._pos = Ipos
    
    @vel.setter
    def vel(self, Ivel):
        self._vel = Ivel

    @property
    def norm_pos(self):
        return np.sqrt(self._pos[-1][0]**2. + self._pos[-1][1]**2.)

    @property
    def norm_vel(self):
        return np.sqrt(self._vel[-1][0]**2. + self._vel[-1][1]**2.)
    
    def append(self, r, v):
        self._pos.append(r)
        self._vel.append(v)


def mplot(p1, p2):
    """matplotlib plot"""
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    axes.set_xlim([-10, 10])
    axes.set_ylim([-10, 10])

    axes.scatter(0, 0, color="red")
    axes.scatter(p1[:, 0], p1[:, 1], color="blue")
    axes.scatter(p2[:, 0], p2[:, 1], color="orange")

    plt.savefig("../assets/orbits.png")


def manimation(p1, p2):
    """matplotlib animation"""
    def _animate(i):
        print(i, len(p1))
        loc_1.set_offsets([p1[i]])
        loc_2.set_offsets([p2[i]])

        if i < 10:
            traj_1.set_offsets(p1[:i+1])
            traj_2.set_offsets(p2[:i+1])
        else:
            traj_1.set_offsets(p1[i-9:i+1])
            traj_2.set_offsets(p2[i-9:i+1])

        return loc_1, loc_2, traj_1, traj_2

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    axes.set_xlim([-10, 10])
    axes.set_ylim([-10, 10])
    loc_0 = axes.scatter([0], [0], label="sun", color="red")
    loc_1 = axes.scatter([], [], label="earth", color="blue")
    loc_2 = axes.scatter([], [], label="jupiter", color="red")
    traj_1 = axes.scatter([], [], s=2, alpha=0.3, color="blue")
    traj_2 = axes.scatter([], [], s=2, alpha=0.3, color="red")

    result  = ani.FuncAnimation(fig, _animate, frames=np.arange(len(p1)), interval=100)
    result.save("../assets/orbits.gif", writer="imagemagick", fps=30, dpi=100)


def simulator():
    # Initialize Planets
    init_earth_pos = [DISTANCE_TO_EARTH/AU_TO_CM, 0]
    init_jupiter_pos= [DISTANCE_TO_JUPITER/AU_TO_CM, 0]
    init_earth_vel = [0, np.sqrt(GG*AU_TO_CM/DISTANCE_TO_EARTH)] 
    init_jupiter_vel = [0, np.sqrt(GG*AU_TO_CM/DISTANCE_TO_JUPITER)]

    earth = Planet(init_earth_pos, init_earth_vel)
    jupiter = Planet(init_jupiter_pos, init_jupiter_vel)

    # CFL timestep
    iterations = 1000 
    h = 0.1/np.max([earth.norm_vel, jupiter.norm_vel])

    for step in range(iterations):
        positions = np.hstack([earth.pos[-1], jupiter.pos[-1]])
        velocities = np.hstack([earth.vel[-1], jupiter.vel[-1]])

        #r, v = integrate_Euler_step(positions, velocities, h)
        #r, v = integrate_EulerCromer_step(positions, velocities, h)
        r, v = integrate_RungeKutta4th_step(positions, velocities, h)

        h = 0.1/np.max([earth.norm_vel, jupiter.norm_vel])

        earth.append(r[:2], v[:2])
        jupiter.append(r[2:], v[2:])

    return earth, jupiter


if __name__=="__main__":
    earth, jupiter = simulator()
    mplot(earth.pos, jupiter.pos)
    manimation(earth.pos, jupiter.pos)
