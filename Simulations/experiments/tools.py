# Code sourced from https://github.com/nruhl25/HorizonCrossings-
# Author: Nathaniel Ruhl

# libraries
import numpy as np
from numpy.linalg import norm
import math as m

# constants
G = 6.6743 * 10 ** (-11)     # Nm^2/kg^2, Gravitational constant
R_EARTH = 6371               # km, known radius of earth
M_EARTH = 5.972 * 10 ** 24   # kg, known mass of Earth
MU = G * M_EARTH             # Nm^2/kg, Earth's gravitational parameter

# Helper functions

def celestial_to_geocentric(alpha, delta):
    x = np.cos(delta)*np.cos(alpha)
    y = np.cos(delta)*np.sin(alpha)
    z = np.sin(delta)
    return np.array([x, y, z])

def period_to_a(T):
    a = (((T ** 2) * G * M_EARTH / (4 * np.pi ** 2)) ** (1. / 3)) / (10 ** 3)  # km
    return a


def a_to_period(a_km):
    a_m = a_km * 10 ** 3   # convert a to meters
    T = np.sqrt((4 * np.pi ** 2 * a_m ** 3) / (G * M_EARTH))   # sec
    return T

def get_Q_matrix(i, raan, aop):
    Q = np.array([[m.cos(raan) * m.cos(aop) - m.sin(raan) * m.sin(aop) * m.cos(i),
                   m.sin(raan) * m.cos(aop) + m.cos(raan) * m.cos(i) * m.sin(aop),
                   m.sin(i) * m.sin(aop)],
                  [-m.cos(raan) * m.sin(aop) - m.sin(raan) * m.cos(i) * m.cos(aop),
                   -m.sin(raan) * m.sin(aop) + m.cos(raan) * m.cos(i) * m.cos(aop),
                   m.sin(i) * m.cos(aop)],
                  [m.sin(raan) * m.sin(i),
                   -m.cos(raan) * m.sin(i),
                   m.cos(i)]])
    return Q

def get_position_from_keplerian(omega, Q, R, T):
    t_array = np.arange(0, T)
    
    perifocal_positions = R * np.array([[np.cos(omega * t) for t in t_array], 
                                    [np.sin(omega * t) for t in t_array], 
                                    [0 for t in t_array]])
    
    ECI_positions = Q.T @ perifocal_positions

    return ECI_positions
