import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
from scipy.signal import find_peaks
import uncertainties as unc
import uncertainties.unumpy as unp

#Funktionen:
def f(x, m, c):
    return m*x+c 

def gauß(x,h,u,s,g):
    return h*np.exp(-((x-u)**2)/(2*s**2))+g


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
plt.figure(constrained_layout=True)
plt.bar(x_eu,N_eu,width=1,label="Messdaten")
plt.plot(peak,N_eu[peak],"rx",label="Peak")
plt.xlabel("Channel")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
#plt.show()


#Energiezuteilung
Energie=[122,245,344,411,444,779,867,964]
x=np.linspace(a,b)

par, cov=curve_fit(f,peak, Energie)
print(cov)
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
c = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
c_f= float(unp.std_devs(par[1]))

print(m,c)

plt.figure(constrained_layout=True)
plt.plot(peak,Energie,"rx",label="Peaks")
plt.plot(x,f(x,m,c),"b-",label="Ausgleichsgerade")
plt.xlabel("Channel")
plt.ylabel("E [keV]")
plt.xlim(a,b)
plt.ylim(0,1000)
plt.legend()
#plt.show()


#Gaußanpassung
d=25

for i in range(len(peak)):

    par, cov=curve_fit(gauß,x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d], p0=[1,peak[i],1,1], sigma=np.sqrt(N_eu[peak[i]-d:peak[i]+d]))
    par = unc.correlated_values(par, cov)
    h = float(unp.nominal_values(par[0]))
    u = float(unp.nominal_values(par[1]))
    s = float(unp.nominal_values(par[2]))
    g = float(unp.nominal_values(par[3]))

    x=np.linspace(x_eu[peak[i]-d],x_eu[peak[i]+d],1000)
    
    plt.figure(constrained_layout=True)
    #plt.errorbar(x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d],yerr=np.sqrt(N_eu[peak[i]-d:peak[i]+d]),fmt="r")
    plt.bar(x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d],width=1,label=f"Messdaten Peak {i+1}")
    plt.plot(x,gauß(x,h,u,s,g),"g-",label="Gauß-Fit")
    plt.legend()
    plt.show()



#------------------------------------------------------------------------------------------

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

