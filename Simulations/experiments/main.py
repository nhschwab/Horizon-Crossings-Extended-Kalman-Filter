# Experimental Script

# libraries
from sim_xray_source import Xray_Source
from hcnm_sim import HCNM_Sim
from LocateR0hc import LocateR0hc

# main function
def main():
    source = Xray_Source(source_name='Example')
    obj = HCNM_Sim(source)
    print(obj.h_unit)
    obj.plot_orbit()

    hc_data = obj.generate_dict()
    r0hc = LocateR0hc(observation_dict=hc_data, earth_shape_string='sphere', r_model_type='circle')
    print(r0hc.r0_hc)
    print(r0hc.r0_2d)
    r0hc.plot_locate_r02d()



if __name__ == "__main__":
    main()