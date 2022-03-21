# Experimental Script

# libraries
import numpy as np
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
    hc_data = obj.generate_dict()

    # locate grazing point of source in orbit
    r0hc = LocateR0hc(observation_dict=hc_data, earth_shape_string='sphere', r_model_type='circle')
    print(f"r0 = {r0hc.r0_hc}")

    # plot orbital trajectory and r0
    obj.plot_orbit(r0hc.r0_hc)

    return None

def crab():

    # simulate crab HC
    source = Xray_Source(source_name='Crab')
    source.ra = np.deg2rad(83.63317) 
    source.dec = np.deg2rad(22.01453)

    # instantiate HCNM Sim object
    obj = HCNM_Sim(source, inclination=51.6, raan=67.99627740036111)
    hc_data = obj.generate_dict()

    # locate grazing point of source in orbit
    r0hc = LocateR0hc(observation_dict=hc_data, earth_shape_string='sphere', r_model_type='circle')
    print(f"r0 = {r0hc.r0_hc}")

    # plot orbital trajectory and r0
    obj.plot_orbit(r0hc.r0_hc)

    return None

def random():

    # simulate random HC
    source = Xray_Source(source_name='Random')

    # instantiate HCNM Sim object
    obj = HCNM_Sim(source)
    hc_data = obj.generate_dict()

    # locate grazing point of source in orbit
    r0hc = LocateR0hc(observation_dict=hc_data, earth_shape_string='sphere', r_model_type='circle')
    print(f"r0 = {r0hc.r0_hc}")

    # plot orbital trajectory and r0
    obj.plot_orbit(r0hc.r0_hc)

    return None





if __name__ == "__main__":
    random()