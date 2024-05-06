import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp

n = np.genfromtxt("EU-Spektrum.txt", unpack=True)

a=0;b=2000 #Intervall
x=np.zeros(len(n[a:b]))
N=np.zeros(len(n[a:b]))

for index, val in enumerate(n[a:b]):
    x[index]=index+a
    N[index]=val


#Plot EU Spektrum
plt.bar(x,N,label="Messdaten")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
plt.show()

