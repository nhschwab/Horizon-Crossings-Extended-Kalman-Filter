# Experimental Simulation of Random X-ray Source

# libraries
import numpy as np 
import matplotlib.pyplot as plt

"""
class: Xray_Source

    generates an object that corresponds to distribution of x-ray photon count rate vs. energy

arguments

    lower_energy: lower energy bound of distribution (default is 0 keV)
    upper_energy: upper energy bound of distribution (default is 10 keV)
    mu: mean of normal distribution (default is 1 keV)
    sigma: standard deviation of distribution (default is 1 keV)
    s: number of samples in distribution (default is 500)


methods

    init: specify spectrum parameters and generate simulated source
    
    show_spectrum: plot histogram of photon energy distribution
    
    show_loc: print the RA and Dec of the source

"""

class Xray_Source():

    def __init__(self, lower_energy=0, upper_energy=10, mu=1, sigma=0.5, s=500):
        
        self.lower_energy = lower_energy
        self.upper_energy = upper_energy
        self.mu = mu
        self.sigma = sigma
        self.s = s

        # upper bound must be greater than lower bound
        if self.upper_energy <= self.lower_energy:
            raise ValueError("Upper energy bound must be greater than lower energy bound")

        # generate spectrum
        self.spectrum = np.abs(np.random.normal(loc=mu, scale=sigma, size=s))
        self.spectrum = np.concatenate((self.spectrum, np.abs(np.random.normal(loc=2, scale=3, size=s))))

        # generate random celestial coordinates
        self.ra = 360 * np.random.ranf()
        self.dec = 180 * np.random.ranf() - 90

    
    def show_spectrum(self):
        plt.hist(self.spectrum, bins=30)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Normalized counts/sec ${keV}^{-1}$")
        plt.title("X-ray Spectrum of Simulated Source")
        plt.show()


    def show_loc(self):
        print(self.ra, self.dec)



if __name__ == "__main__":
    instance = Xray_Source()
    instance.show_loc()
    instance.show_spectrum()