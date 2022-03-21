# Simulation of Orbital Model and Horizon Crossing

# libraries
import numpy as np 
import matplotlib.pyplot as plt
import tools

from sim_xray_source import Xray_Source

'''
class: HCNM_Sim

    generates a simulated oribital model with artifical horizon crossings

    User instantiates an object of the class by passing an X-ray source, which is an object o the class Xray_Source

    The user can manually alter the attributes of object or leave as is for random orbital parameters

    Once satisfied with orbital parameters, user can execute generate_dict command, which will produce
    the python dict object that can be inputted to the HCNM method for trial and analysis

arguments:

    sources: array Xray_Source objects
    
'''

class HCNM_Sim():

    def __init__(self, source):

        self.source_name = source.source_name

        # ORBITAL PARAMETERS
        
        # orbital radius
        self.R_orbit = 6803.76 # km
        # orbital period
        self.T = tools.a_to_period(self.R_orbit)
        # angular velocity
        self.OMEGA_ORB = 2 * np.pi / self.T

        # KEPLERIAN ELEMENTS
        self.inclination = 51.6 # 180 * np.random.ranf()
        self.raan = 360 * np.random.ranf()
        self.aop = 0 # CIRCULAR ORBIT
        
        # the z-component of the pole vector is defined as the cos of the inclination
        # the other components are arbitrarily defined
        a = np.random.ranf()
        b = np.sqrt(1 - np.cos(self.inclination)**2 - a**2)
        self.h_unit = np.array([a, b, np.cos(self.inclination)])

        # the x-component of the line of nodes is defined as the cos of the raan
        # the other components are arbitrarily defined
        self.n_unit = np.array([np.cos(self.raan), np.random.ranf(), np.random.ranf()])

        # Q matrix: transformation from perifical frame to planet-centered coordinate frame
        self.Q = tools.get_Q_matrix(self.inclination, self.raan, self.aop)

        self.positions = tools.get_position_from_keplerian(self.OMEGA_ORB, self.Q, self.R_orbit, self.T)

        # information about the source necessary for geometric computation
        self.starECI = tools.celestial_to_geocentric(source.ra, source.dec)
        self.starECI_proj = tools.proj_on_orbit(self.starECI, self.h_unit)


    # method to plot orbital model
    def plot_orbit(self):
    
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.positions[0], self.positions[1], self.positions[2])
        plt.show()

        return None

    # method generates a dict object that is the necessary input for LocateR0hc
    # this method should only be called once all orbital parameters are set for the simulation
    def generate_dict(self):
        d = {'Source_Name': self.source_name}
        
        d['T'] = self.T
        d['h_unit'] = self.h_unit
        d['R_orbit'] = self.R_orbit
        d['OMEGA_ORB'] = self.OMEGA_ORB
        d['Q'] = self.Q
        d['starECI_proj'] = self.starECI_proj
        d['starECI'] = self.starECI

        return d



if __name__ == "__main__":
    source = Xray_Source('simulated source')
    obj = HCNM_Sim(source)
    obj.plot_orbit()





