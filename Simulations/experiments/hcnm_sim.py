# Simulation of Orbital Model and Horizon Crossing

# libraries
import numpy as np 
import matplotlib.pyplot as plt
import tools
from orbital import earth, KeplerianElements, plot3d, earth_sidereal_day
from scipy.constants import kilo

from sim_xray_source import Xray_Source

'''
class: HCNM_Sim

    generates a simulated oribital model with artifical horizon crossings

    User instantiates an object of the class by passing an array of sources

    The user can manually alter the attributes of object or leave as is for random orbital parameters

    Once satisfied with orbital parameters, user can execute generate_dict command, which will produce
    the python dict object that can be inputted to the HCNM method for trial and analysis

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
        self.T = tools.a_to_period(self.R_orbit)
        # angular velocity
        self.OMEGA_ORG = 2 * np.pi / self.T

        # KEPLERIAN ELEMENTS
        self.inclination = 180 * np.random.ranf()
        self.raan = 360 * np.random.ranf()
        self.aop = 0 # CIRCULAR ORBIT
        
        # the z-component of the pole vector is defined as the cos of the inclination
        # the other components are arbitrarily defined
        h = np.array([np.random.ranf(), np.random.ranf(), np.cos(self.inclination)])
        self.h_unit = h / np.sqrt(np.sum(h**2)) # normalize the vector

        # the x-component of the line of nodes is defined as the cos of the raan
        # the other components are arbitrarily defined
        self.n_unit = np.array([np.cos(self.raan), np.random.ranf(), np.random.ranf()])

        # Q matrix: transformation from perifical frame to planet-centered coordinate frame
        self.Q = tools.get_Q_matrix(self.inclination, self.raan, self.aop)

        self.positions = tools.get_position_from_keplerian(self.OMEGA_ORG, self.Q, self.R_orbit, self.T)


    # method to plot orbital model
    def plot_orbit(self):
    
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.positions[0], self.positions[1], self.positions[2])
        plt.show()


if __name__ == "__main__":
    obj = HCNM_Sim([])
    print(obj.plot_orbit())





