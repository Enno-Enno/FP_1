import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp

#Aus Plot1.py:
m=0.20734414185681416

def WQ_compton(alpha,theta):
    r_0, _, delta_r = const["classical electron radius"]
    Z = 40
    cos = lambda x: np.cos(x)
    WQ = (
        Z
        * r_0**2
        * (1 / (alpha * (1 - cos(theta)) + 1)) ** 2
        * (cos(theta) ** 2 + 1) / 2
        * (
            (alpha**2 * (1 - cos(theta)) ** 2)
            / (((alpha * (-cos(theta) + 1) + 1) * (cos(theta) ** 2 + 1)))
            + 1
        )
    )  # - Derivative(sigma, Omega) == 0
    return WQ


n_cs = np.genfromtxt("137Cs-Spektrum.txt", unpack=True)

a=0;b=4000 #Intervall
x_cs=np.zeros(len(n_cs[a:b]))
N_cs=np.zeros(len(n_cs[a:b]))

for index, val in enumerate(n_cs[a:b]):
    x_cs[index]=(index+a)*m
    N_cs[index]=val

plt.figure(constrained_layout=True)
plt.bar(x_cs,N_cs,width=1,label="Messdaten 137Cs")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
#plt.show()


#------------------------------------------------------------------------------------------------

n_ba = np.genfromtxt("133Ba-Spektrum.txt", unpack=True)

a=0;b=3000 #Intervall
x_ba=np.zeros(len(n_ba[a:b]))
N_ba=np.zeros(len(n_ba[a:b]))

for index, val in enumerate(n_ba[a:b]):
    x_ba[index]=(index+a)*m
    N_ba[index]=val

plt.figure(constrained_layout=True)
plt.bar(x_ba,N_ba,width=1,label="Messdaten 133Ba")
plt.xlabel("Energie E [kev]")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
#plt.show()


#------------------------------------------------------------------------------------------------

n_un = np.genfromtxt("000Unbekannt-Spektrum.txt", unpack=True)

a=0;b=8000 #Intervall
x_un=np.zeros(len(n_un[a:b]))
N_un=np.zeros(len(n_un[a:b]))

for index, val in enumerate(n_un[a:b]):
    x_un[index]=(index+a)*m
    N_un[index]=val

plt.figure(constrained_layout=True)
plt.bar(x_un,N_un,width=1,label="Messdaten Unbekannt")
plt.xlabel("Energie E [kev]")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
#plt.show()
