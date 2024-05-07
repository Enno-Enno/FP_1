import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
from scipy.signal import find_peaks
import uncertainties as unc
import uncertainties.unumpy as unp

#Funktionen:
def g(x, m, c):
    return m*x+c 



n_eu = np.genfromtxt("152Eu-Spektrum.txt", unpack=True)

a=0;b=5000 #Intervall
x_eu=np.zeros(len(n_eu[a:b]))
N_eu=np.zeros(len(n_eu[a:b]))

for index, val in enumerate(n_eu[a:b]):
    x_eu[index]=index+a
    N_eu[index]=val

peak,_=find_peaks(N_eu,height=100)
peak2,_=find_peaks(N_eu[int(b/3):b],height=25,distance=10)
peak2=peak2+int(b/3)
peak=np.concatenate([peak, peak2])
print(peak)
print(N_eu[peak])

#Plot EU Spektrum
plt.bar(x_eu,N_eu,width=1,label="Messdaten")
plt.plot(peak,N_eu[peak],"rx",label="Peak")
#plt.plot(peak2,N[peak2],"rx")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
plt.show()

Energie=[122,245,344,411,444,779,867,964]

x=np.linspace(a,b)

par, cov=curve_fit(g,peak, Energie)
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
c = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
c_f= float(unp.std_devs(par[1]))

print(m,c)

plt.plot(peak,Energie,"rx",label="Peaks")
plt.plot(x,g(x,m,c),"b-",label="Ausgleichsgerade")
plt.xlabel("Channel")
plt.ylabel("E [keV]")
plt.xlim(a,b)
plt.ylim(0,1000)
plt.legend()
plt.show()





#------------------------------------------------------------------------------------------

n_cs = np.genfromtxt("137Cs-Spektrum.txt", unpack=True)

a=0;b=4000 #Intervall
x_cs=np.zeros(len(n_cs[a:b]))
N_cs=np.zeros(len(n_cs[a:b]))

for index, val in enumerate(n_cs[a:b]):
    x_cs[index]=index+a
    N_cs[index]=val

plt.bar(x_cs,N_cs,width=1,label="Messdaten")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
#plt.show()


#------------------------------------------------------------------------------------------------

n_ba = np.genfromtxt("133Ba-Spektrum.txt", unpack=True)

a=0;b=3000 #Intervall
x_ba=np.zeros(len(n_ba[a:b]))
N_ba=np.zeros(len(n_ba[a:b]))

for index, val in enumerate(n_ba[a:b]):
    x_ba[index]=index+a
    N_ba[index]=val

plt.bar(x_ba,N_ba,width=1,label="Messdaten")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
#plt.show()


#------------------------------------------------------------------------------------------------

n_un = np.genfromtxt("000Unbekannt-Spektrum.txt", unpack=True)

a=0;b=8000 #Intervall
x_un=np.zeros(len(n_un[a:b]))
N_un=np.zeros(len(n_un[a:b]))

for index, val in enumerate(n_un[a:b]):
    x_un[index]=index+a
    N_un[index]=val

plt.bar(x_un,N_un,width=1,label="Messdaten")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
#plt.show()

