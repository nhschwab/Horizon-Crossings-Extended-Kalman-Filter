# Experimental Script

# libraries
from sim_xray_source import Xray_Source

# main function
def main():
    instance = Xray_Source()
    instance.show_loc()
    instance.show_spectrum()


if __name__ == "__main__":
    main()