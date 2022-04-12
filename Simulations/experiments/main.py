# Experimental Script

# libraries
import numpy as np
from astropy.io import ascii
from sim_xray_source import Xray_Source
from hcnm_sim import HCNM_Sim
from LocateR0hc import LocateR0hc

# simulate v4641
def v4641():

    # simulate V4641 Feb 3 HC
    source = Xray_Source(source_name='V4641_Feb3')
    source.ra = np.deg2rad(274.839) 
    source.dec = np.deg2rad(-25.407)

    # instantiate HCNM Sim object
    obj = HCNM_Sim(source, inclination=51.6, raan=293.1)

    # generate simulated event file
    data = obj.generate_EVT()

    # write event file
    ascii.write(data, 'v4641_events.dat')

    return None

def crab():

    # simulate crab HC
    source = Xray_Source(source_name='Crab')
    source.ra = np.deg2rad(83.63317) 
    source.dec = np.deg2rad(22.01453)

    # instantiate HCNM Sim object
    obj = HCNM_Sim(source, inclination=51.7, raan=68.0)

    # generate simulated event file
    data = obj.generate_EVT()

    # write event file
    ascii.write(data, 'crab_events.dat')

    return None

def random():

    # simulate random HC
    source = Xray_Source(source_name='Random')

    # instantiate HCNM Sim object
    obj = HCNM_Sim(source)

    # generate simulated event file
    data = obj.generate_EVT()

    # write event file
    ascii.write(data, 'sim_events3.dat')

    return None


if __name__ == "__main__":
    v4641()