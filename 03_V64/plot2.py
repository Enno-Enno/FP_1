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
    return  4*np.pi/(lamVac) *T_0/(np.cos(theta_rad))* theta_rad**2

def Rcombined(theta_deg):
    theta_rad = degReturnRad(theta_deg)
    theta_01_deg = 10 
    theta_02_deg = -10
    theta_01_rad = degReturnRad(theta_01_deg) 
    theta_02_rad = degReturnRad(theta_02_deg) 
    return (R(theta_01_rad + theta_rad) + R(theta_rad + theta_02_rad))  


def nGlass(DeltaPhi):
    return 1/(1-(2*DeltaPhi)/(Rcombined(10)- Rcombined(0))) 

def getk(n,theta):
    DeltaPhi =  np.abs(Rcombined(theta)- Rcombined(0)) * (n-1)/(2*n)
    return DeltaPhi /(2*np.pi)

# for val in counter_glass:
    # print(val, r"\\ ")

uCounterGlass = unc.ufloat(np.mean(counter_glass), np.std(counter_glass))
print(uCounterGlass)
DeltaPhi = 2* np.pi* uCounterGlass
n = nGlass(DeltaPhi)
print(n)
print("Kontrolle", getk(n, 10))

### Refractive Index air
print(d.counter_air)
difference_array = np.zeros((np.shape(d.counter_air)[0],np.shape(d.counter_air)[1]-1))
for i,_ in enumerate(difference_array[0,:]):
    difference_array[0,i] = d.counter_air[0,i+1] - d.counter_air[0,i]
    difference_array[1,i] = d.counter_air[1,i+1] - d.counter_air[1,i]
    difference_array[2,i] = d.counter_air[2,i+1] - d.counter_air[2,i]
    difference_array[3,i] = d.counter_air[3,i+1] - d.counter_air[3,i]

dcounter = unc.ufloat(np.mean(difference_array),np.std(difference_array))
print(dcounter)
endarray = d.counter_air[:,-1]
end = unc.ufloat(np.mean(endarray),np.std(endarray))
print(end)

def deltaN(DeltaPhi_rad, L):
    return DeltaPhi_rad/(2*np.pi*L) * lamVac
L = unc.ufloat(100, 0.1) *1e-3 #m
delN = deltaN(end*2*np.pi, L)
print("Delta n ", delN)
delNFine = deltaN(dcounter*np.pi, L) * len(difference_array[0,:])
print(delNFine)

for i,_ in enumerate(d.counter_air[0,:]):
    print(d.counter_air[0,i], " & ",d.counter_air[1,i], " & ",d.counter_air[2,i], " & ",d.counter_air[3,i], r" \\",)