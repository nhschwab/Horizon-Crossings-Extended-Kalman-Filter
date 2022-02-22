# Simulation of Orbital Model and Horizon Crossing

# libraries
import numpy as np 
import matplotlib.pyplot as plt
import tools

from sim_xray_source import Xray_Source

'''
class: HCNM_Sim

    generates a simulated oribital model with artifical horizon crossings

arguments:

    sources: array Xray_Source objects
    
'''

class HCNM_Sim():

    def __init__(self, sources):

        self.sources = sources
        
        # number of considered sources
        self.N = len(self.sources)

        # ORBITAL PARAMETERS
        
        # orbital radius
        self.R_orbit = 6803.76 # km
        # orbital period
        self.T = tools.a_to_period(self.R)
        # angular velocity
        self.OMEGA_ORG = 2 * np.pi / self.T

        # KEPLERIAN ELEMENTS
        self.inclination = 180 * np.random.randf()
        self.raan = 360 * np.random.ranf()
        self.aop = 360 * np.random.ranf()

