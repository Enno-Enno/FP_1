import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import data as d
import pprint
from plot1 import degReturnRad
import numpy.polynomial.polynomial as pol

pp = pprint.PrettyPrinter(indent=4)
lamVac = 632.990 * 1e-9# m
counter_glass = d.RefIndexGlasNumbers
T_0 = 1e-3 #m

# Falsche Formeln
# def T(theta_rad):
    # '''return thickness of glass for the laser'''
    # return T_0/np.cos(theta_rad)
    # return T_0
# 
    # 
# def R(theta_rad):
    # return  2*np.pi/(lamVac) *T(theta_rad) * theta_rad**2
# 
# bis hier
def Rcombined(theta_deg):
    theta_rad = degReturnRad(theta_deg)
    theta_0_deg = 10 
    # theta_02_deg = -10
    theta_0_rad = degReturnRad(theta_0_deg) 
    # theta_02_rad = degReturnRad(theta_02_deg) 
    # return (R(theta_01_rad + theta_rad) + R(theta_rad + theta_02_rad))  
    return (4*np.pi *T_0)/(lamVac) * theta_rad * theta_0_rad




def nGlass(DeltaPhi):
    return 1/(1-(DeltaPhi)/(Rcombined(10)))# - Rcombined(0))) 

def getk(n,theta):
    DeltaPhi =  np.abs(Rcombined(theta)- Rcombined(0)) * (n-1)/(2*n)
    return DeltaPhi /(np.pi)

# for val in counter_glass:
    # print(val, r"\\ ")

uCounterGlass = unc.ufloat(np.mean(counter_glass), np.std(counter_glass))
print(uCounterGlass)
DeltaPhi = 2 *np.pi* uCounterGlass
n = nGlass(DeltaPhi)
print(n)
print("Kontrolle", getk(n, 10))

print("Refractive Index air ----------------------------------")
print(d.counter_air)
# difference_array = np.zeros((np.shape(d.counter_air)[0],np.shape(d.counter_air)[1]-1))
# for i,_ in enumerate(difference_array[0,:]):
    # difference_array[0,i] = d.counter_air[0,i+1] - d.counter_air[0,i]
    # difference_array[1,i] = d.counter_air[1,i+1] - d.counter_air[1,i]
    # difference_array[2,i] = d.counter_air[2,i+1] - d.counter_air[2,i]
    # difference_array[3,i] = d.counter_air[3,i+1] - d.counter_air[3,i]
# 
# dcounter = unc.ufloat(np.mean(difference_array),np.std(difference_array))
# print(dcounter)
endarray = d.counter_air[:,-1]
end = unc.ufloat(np.mean(endarray),np.std(endarray))
print(end)
# 
def deltaN(DeltaPhi_rad, L):
    return DeltaPhi_rad/(2*np.pi*L) * lamVac
L = unc.ufloat(100, 0.1) *1e-3 #m
delN = deltaN(end*2*np.pi, L)
print("Delta n ", delN)
# delNFine = deltaN(dcounter*np.pi, L) * len(difference_array[0,:])
# print(delNFine)
# 
for i,_ in enumerate(d.counter_air[0,:]):
    print(d.pressures_mbar[i], " & ",d.counter_air[0,i], " & ",d.counter_air[1,i], " & ",d.counter_air[2,i], " & ",d.counter_air[3,i], r" \\",)
# 
# delN = deltaN(d.counter_air.flatten()*2*np.pi,L)
# nomN = unp.nominal_values(delN)
# 
# print(np.shape(nomN))
deltaPhis= d.counter_air.flatten()*2*np.pi
pressures_mbar = np.append(d.pressures_mbar,[d.pressures_mbar,d.pressures_mbar,d.pressures_mbar])
par, cov = np.polyfit(pressures_mbar,unp.nominal_values(deltaN(deltaPhis,L)),1, cov=True)

par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
d1 = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
d1_f= float(unp.std_devs(par[1]))

A_air = None
alphamol = unc.ufloat(2.118, 0.091)*e-29 # m³

n_lorentz = lambda T, p:np.sqrt((1 + (A_air*P)/(R*T))/(1 -(A_air*P)/(R*T))) 


print(par[0], par[1])
xplot= np.linspace(0,1000)

plt.figure(constrained_layout=True)
plt.plot(pressures_mbar, unp.nominal_values(deltaN(deltaPhis,L)),".",label="measured values")
plt.plot(xplot, m*xplot+ d1, label = "linear fit")

plt.xlabel(r"$p/\text{mBar}$")
plt.ylabel(r" $\Delta n = n-1$")
plt.legend()
plt.savefig("build/plot2.pdf")
# plt.show()