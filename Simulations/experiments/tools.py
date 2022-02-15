# Code sourced from https://github.com/nruhl25/HorizonCrossings-
# Author: Nathaniel Ruhl

# libraries
import numpy as np
from numpy.linalg import norm
from openpyxl import Workbook
from openpyxl import load_workbook
import os.path
import datetime
import pymap3d as pm
import numbers

# Helper functions

def celestial_to_geocentric(alpha, delta):
    x = np.cos(delta)*np.cos(alpha)
    y = np.cos(delta)*np.sin(alpha)
    z = np.sin(delta)
    return np.array([x, y, z])
