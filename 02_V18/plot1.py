import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
from scipy.signal import find_peaks
import uncertainties as unc
import uncertainties.unumpy as unp

n = np.genfromtxt("EU-Spektrum.txt", unpack=True)

a=0;b=3000 #Intervall
x=np.zeros(len(n[a:b]))
N=np.zeros(len(n[a:b]))

for index, val in enumerate(n[a:b]):
    x[index]=index+a
    N[index]=val

peak,_=find_peaks(N,height=100)
print(peak)
print(N[peak])

#Plot EU Spektrum
plt.bar(x,N,width=1,label="Messdaten")
plt.plot(peak,N[peak],"rx",label="Peak")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
plt.show()

