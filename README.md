This repository contains a library for simulating horizon crossing events as detected by an X-ray telescope aboard a satellite.

Author: Noah Schwab

Advisor: Andrea Lommen

Haverford College, Department of Physics and Astronomy

How to use this repository:

To start, clone this repository to local using standard git protocols. 

To run a simulation of a given satellite's orbit, run main.py from the command line. Within the main.py script, you will see various functions. These functions write a new .dat file containing the artificial X-ray photon detection data. The .dat files are a data table file extension that has a form similar to the .mkf and .evt files provided by HEASARC. The column "TIME" corresponds to seconds since the start of the simulation, and the column "PI" corresponds to energy of the X-ray photon detected at that time stamp. A value of 0 in the "PI" column signifies that no detection was made at that time.

In order to select parameters for a horizon crossing simulation, write a new function in main.py given thetemplate provided by existing functions. This template enables the user to specificy the location of the simulated X-ray source, as well as the orbital parameters of the simulated satellite. This is accomplished via the classes HCNM_Sim and Xray_Source, which are written in the python scripts hcnm_sim.py and sim_xray_source.py, respectively. 

Once you have written a function in main.py according to the template with the appropriate parameters for a simulated horizon crossing of interest, you can make that function executable by writing beneath the line 'if __name__ == "__main__":'. Now when you run main.py from the command line, the program will output a .dat file with the simulated X-ray photon detection data in your current directory. 

In order to analyze the simulated horizon crossing data, you must use Nathaniel Ruhl's repository, "HorizonCrossings-". Several of the scripts and tools in Nathaniel's repo can be used to analyze the simulated data, however the most helpful for visualizing the horizon crossing is Python_Scripts/HCNM/TransmitOrbitModel.py. To use this script with the simulated data, you must create a python scripts that includes a dictionary with predefined key, value pairs, according to the template of crab03_dict.py and feb3_dict.py found in the Observations directory under their respective subdirectories. In these dictionaries you can read in the simulated X-ray photon detections according to the protocols dictated by Nathaniel's repo. By modifying the scripts that are read-in and defined in TransmitOrbitModel.py, as well as determining the relevant functions to include in the executable of the script, you can run TransmitOrbitModel.py from the command line and see different visualizations of the simulated horizon crossing.

Notes:

The parameters passed to the simulator must be carefully chosen such that the geometry is possible for detecting a horizon crossing. If random parameters are passed for the location of the X-ray source and the orbital parameters, there is a high probability that an error will be raised due to geometric incompatabilities, or the simulated X-ray photon detection data will not correspond to a horizon crossing. 

Also note that in many places in the code, files are being read in and specified paths correspond to my own local directory structure. You can remedy this by simply changing the string passed in those lines to your own local path names.




