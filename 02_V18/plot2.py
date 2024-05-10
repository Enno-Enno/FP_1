import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks


#Aus Plot1.py:
m=0.20734414185681416

#Funktionen:

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

def gauß(x,h,u,s,g):
    return h*np.exp(-((x-u)**2)/(2*s**2))+g



n_cs = np.genfromtxt("137Cs-Spektrum.txt", unpack=True)

a=0;b=4000 #Intervall
x_cs=np.zeros(len(n_cs[a:b]))
N_cs=np.zeros(len(n_cs[a:b]))

for index, val in enumerate(n_cs[a:b]):
    x_cs[index]=(index+a)*m
    N_cs[index]=val
    
    
peak,_=find_peaks(N_cs,height=70)
peak2,_=find_peaks(N_cs[2200:2500],height=42,distance=10)
peak2=peak2+2200
peak=np.concatenate([peak, peak2])
print(m*peak)
print(N_cs[peak])


plt.figure(constrained_layout=True)
plt.bar(x_cs,N_cs,width=m,label="Messdaten 137Cs")
plt.plot(m*peak,N_cs[peak],"rx",label="Peak")
plt.xlabel("Energie")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
plt.show()

#Zentelbreite und Halbwertsbreite
d=15

par, cov=curve_fit(gauß,x_cs[peak[4]-d:peak[4]+d],N_cs[peak[4]-d:peak[4]+d], p0=[12,m*peak[4],1,1] ,sigma=np.sqrt(N_cs[peak[4]-d:peak[4]+d])+0.01)
par = unc.correlated_values(par, cov)
h = float(unp.nominal_values(par[0]))
u = float(unp.nominal_values(par[1]))
s = float(unp.nominal_values(par[2]))
g = float(unp.nominal_values(par[3]))
h_f= float(unp.std_devs(par[0]))
u_f= float(unp.std_devs(par[1]))
s_f= float(unp.std_devs(par[2]))
g_f= float(unp.std_devs(par[3]))

h1=unp.uarray(h,h_f)
u1=unp.uarray(u,u_f)
s1=unp.uarray(abs(s),s_f)
g1=unp.uarray(g,g_f)

delta_xh=s1*unp.sqrt(s1*unp.log(2*h1/(h1-g1)))

x=np.linspace(x_cs[peak[4]-d],x_cs[peak[4]+d],1000)

plt.figure(constrained_layout=True)
plt.bar(x_cs[peak[4]-d:peak[4]+d],N_cs[peak[4]-d:peak[4]+d],width=m,yerr=np.sqrt(N_cs[peak[4]-d:peak[4]+d]),label=f"Messdaten Photopeak ")
plt.plot(x,gauß(x,h,u,s,g),"g-",label="Gauß-Fit")
plt.plot(x,gauß(u+delta_xh,h,u,s,g),"b--")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.xlim(x_cs[peak[4]-d],x_cs[peak[4]+d])
plt.legend()
plt.show()





#------------------------------------------------------------------------------------------------

n_ba = np.genfromtxt("133Ba-Spektrum.txt", unpack=True)

a=0;b=3000 #Intervall
x_ba=np.zeros(len(n_ba[a:b]))
N_ba=np.zeros(len(n_ba[a:b]))

for index, val in enumerate(n_ba[a:b]):
    x_ba[index]=(index+a)*m
    N_ba[index]=val

plt.figure(constrained_layout=True)
plt.bar(x_ba,N_ba,width=m,label="Messdaten 133Ba")
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
plt.bar(x_un,N_un,width=m,label="Messdaten Unbekannt")
plt.xlabel("Energie E [kev]")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
#plt.show()
