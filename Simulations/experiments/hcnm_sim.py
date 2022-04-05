# Simulation of Orbital Model and Horizon Crossing

# libraries
import numpy as np 
import matplotlib.pyplot as plt
import tools
from astropy.table import Table
from scipy.interpolate import interp1d as interp1d

from sim_xray_source import Xray_Source
from LocateR0hc import LocateR0hc

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

    def __init__(self, source, inclination=None, raan=None):

        self.source_name = source.source_name

        # ORBITAL PARAMETERS
        
        # orbital radius
        self.R_orbit = 6803.76 # km
        # orbital period
        self.T = tools.a_to_period(self.R_orbit)
        # angular velocity
        self.OMEGA_ORB = 2 * np.pi / self.T

        # KEPLERIAN ELEMENTS
        if inclination==None:
            self.inclination =  180 * np.random.ranf()
        else:
            self.inclination = np.deg2rad(inclination)

        if raan==None:
            self.raan = 360 * np.random.ranf()
        else:
            self.raan = np.deg2rad(raan)

        self.aop = 0 # CIRCULAR ORBIT

        # Q matrix: transformation from perifical frame to planet-centered coordinate frame
        self.Q = tools.get_Q_matrix(self.inclination, self.raan, self.aop)

        # normalized pole vector is the third row of the Q matrix
        self.h_unit = self.Q[2]

        self.positions = tools.get_position_from_keplerian(self.OMEGA_ORB, self.Q, self.R_orbit, self.T)

        # information about the source necessary for geometric computation
        self.starECI = tools.celestial_to_geocentric(source.ra, source.dec)
        self.starECI_proj = tools.proj_on_orbit(self.starECI, self.h_unit)


    # method to plot orbital model and grazing point of source
    def plot_orbit(self, r0):
    
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.positions[0], self.positions[1], self.positions[2])
        ax.scatter(r0[0], r0[1], r0[2], c='red', label='r0')
        ax.set_title(f"Simulated Orbital Trajectory and Computed r0 of {self.source_name}")
        ax.set_xlabel("km")
        ax.set_ylabel("km")
        ax.set_zlabel("km")
        plt.legend()
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

    # method to generate a table consisting of times, positions, velocities 
    # of simualted orbital trajectory
    # output serves as simulated MKF file
    def generate_MKF(self):
        data = Table()
        data['TIME'] = np.arange(0, self.T)
        data['POSITION'] = self.positions.T
        data['VELOCITY'] = (self.positions * self.OMEGA_ORB).T
        return data

    
    # method to generate a table consisting of times and x-ray photon observation energy values
    # output serves as simulated EVT file
    def generate_EVT(self):

        # read in atmospheric densities from msis model
        msis_data = np.array([x.split(',') for x in open('/Users/noahhschwab/Desktop/Thesis Work/Horizon-Crossings-Extended-Kalman-Filter/Simulations/experiments/msis00_v4641_feb3.txt').readlines()])
        msis_alt_list = msis_data[:, 0].astype('float64')
        msis_density_list = msis_data[:, 1].astype('float64')

        # use scipy interp1d to get function for evaluating at a given altitude
        f = interp1d(msis_alt_list, msis_density_list, kind='cubic')

        hc_data = self.generate_dict()
        
        # time array, define as t0 + 300 seconds
        r0hc = LocateR0hc(observation_dict=hc_data, earth_shape_string='sphere', r_model_type='circle')
        ind = np.argmin(np.sum((self.positions.T - r0hc.r0_hc)**2, axis=1))
        t_array = np.arange(0, self.T)[ind: ind+300]

        # instantiate empty pi array
        pi_array = []

        # step size in los
        ds = 0.5 # km

        # at each time step, generate a photon taken randomly from x-ray spectrum of simulated source
        # FOR NOW, we are treating x-ray spectrum as gaussian distribution from 0 keV to 5 keV
        for i, t in enumerate(t_array):

            photon = np.random.normal(loc=2.5, scale=1)

            # count photon as an observation of its randomly assigned value [0, 1] falls
            # below the optical depth value according to Beer's Law
            rand_val = np.random.rand()

            # compute LOS vector at each time
            los = tools.line_of_sight(self.positions.T[ind: ind+300][i], self.starECI, ds)
            
            # compute radial altitude of point on LOS vector 
            los_mag = np.sqrt(los[:, 0]**2 + los[:, 1]**2 + los[:, 2]**2)
            altitude_list = los_mag - np.ones_like(los_mag) * tools.R_EARTH
            
            # last index of integration is when altitude is minimum, corresponds to half LOS
            end_ind = np.argmin(altitude_list)

            # redefine altitude list associated with half los
            altitude_list = altitude_list[:end_ind]

            # compute atmopsheric densities at each altitude in altitude_list from msis model
            density_list = f(altitude_list)

            # compute absorption probability from Beer's Law using numerical integration
            optical_depth = np.sum([density_list * cross_section * ds])
            tau = 2 * optical_depth
            print(tau)

            
            
            
            





if __name__ == "__main__":
    source = Xray_Source('simulated source')
    obj = HCNM_Sim(source)
    print(obj.generate_MKF())
    obj.generate_EVT()





