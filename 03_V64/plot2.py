import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import data as d
import pprint
from plot1 import degReturnRad

pp = pprint.PrettyPrinter(indent=4)
lamVac = 632.990 * 1e-9# m
counter_glass = d.RefIndexGlasNumbers
T_0 = 1e-3 #m

def T(theta_rad):
    '''return thickness of glass for the laser'''
    return T_0/np.cos(theta_rad)
    ...

    
def R(theta_rad):
    return  2*np.pi/(lamVac) *T_0/(np.cos(theta_rad))* theta_rad**2

def Rcombined(theta_deg):
    theta_rad = degReturnRad(theta_deg)
    theta_01_deg = 10 
    theta_02_deg = -10
    theta_01_rad = degReturnRad(theta_01_deg) 
    theta_02_rad = degReturnRad(theta_02_deg) 
    return (R(theta_01_rad+ theta_rad) + R(theta_rad + theta_02_rad))  


def nGlass(DeltaPhi):
    return 1/(1-(2*DeltaPhi)/(Rcombined(10)- Rcombined(0))) 

# for val in counter_glass:
    # print(val, r"\\ ")

uCounterGlass = unc.ufloat(np.mean(counter_glass), np.std(counter_glass))
print(uCounterGlass)
DeltaPhi = 2*np.pi* uCounterGlass
n = nGlass(DeltaPhi)
print(n)
