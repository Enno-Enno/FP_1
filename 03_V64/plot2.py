import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import data as d
import pprint
from plot1 import degReturnRad

pp = pprint.PrettyPrinter(indent=4)

def nGlass(k, rotation_angle_deg):
    # TODO Richtige Formel Erstellen mit Verkippung des Glases
    lamVac = 632.990 # nm
    theta = degReturnRad(rotation_angle_deg)
    return 1/(1-k* lamVac* np.cos(theta)/(theta**2))#  k = Delta Phi / (2pi)

counter_glass = d.RefIndexGlasNumbers
