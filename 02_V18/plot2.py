import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
import fktn ## Enthält die Hilfreichen Funktionen fürs Compton Fitten :)


#Aus Plot1.py:
m=0.20734414185681416

#Funktionen:
## In anderer Datei!


WQ_compton = lambda E_gamma,theta: fktn.WQ_compton

def gauß(x,h,u,s,g):
    return h*np.exp(-((x-u)**2)/(2*s**2))+g

def f(x, m, c):
    return m*x+c 


n_cs = np.genfromtxt("137Cs-Spektrum.txt", unpack=True)

a=0;b=4000 #Intervall
x_cs=np.zeros(len(n_cs[a:b]))
N_cs=np.zeros(len(n_cs[a:b]))

for index, val in enumerate(n_cs[a:b]):
    x_cs[index]=(index+a)*m ## Interpretiere ich als Energie. Einheit keV?
    N_cs[index]=val
    
    
peak,_=find_peaks(N_cs,height=70)
peak2,_=find_peaks(N_cs[2200:2500],height=42,distance=10)
peak2=peak2+2200
peak=np.concatenate([peak, peak2])
print("Peak Energien Cs",m*peak)
print(N_cs[peak])


plt.figure(constrained_layout=True)
plt.bar(x_cs,N_cs,width=m,label="Messdaten 137Cs")
plt.plot(fktn.WQ_Energie(x_cs*20,662.46453323,30, 10),"r", label= "Compton \"fit\"")
plt.plot(m*peak,N_cs[peak],"rx",label="Peak")
plt.xlabel("Energie")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
plt.show()

#Zentelbreite und Halbwertsbreite Inhalt Gauß
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

delta_xH=s*np.sqrt(2*np.log(2*h/(h-g))) #Halbwertsbreite
xH=np.linspace(u-delta_xH,u+delta_xH,1000)

delta_xZ=s*np.sqrt(2*np.log(10*h/(h-g))) #Zentelbreite
xZ=np.linspace(u-delta_xZ,u+delta_xZ,1000)


x=np.linspace(x_cs[peak[4]-d],x_cs[peak[4]+d],1000)

plt.figure(constrained_layout=True)
plt.bar(x_cs[peak[4]-d:peak[4]+d],N_cs[peak[4]-d:peak[4]+d],width=m,yerr=np.sqrt(N_cs[peak[4]-d:peak[4]+d]),label=f"Messdaten Photopeak ")
plt.plot(x,gauß(x,h,u,s,g),"g-",label="Gauß-Fit")
plt.plot(xH,f(xH,0.000001,0.5*h),"y--", label="Halb Peak Breite")
plt.plot(xZ,f(xZ,0.000001,0.1*h),"r--", label="Zehntel Peak Breite")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.xlim(u-3,u+3)
plt.legend()
plt.show()

h=unp.uarray(h,h_f)
u=unp.uarray(u,u_f)
s=unp.uarray(abs(s),s_f)
g=unp.uarray(g,g_f)

print("Zentelbreite=",2*s*unp.sqrt(2*unp.log(10*h/(h-g))))
print("Halbwertsbreite=",2*s*unp.sqrt(2*unp.log(2*h/(h-g))))
print("Inhalt=",np.sqrt(2*np.pi)*h*s)

print("\n",h, "\n ", u ,"\n" ,s ,"\n", g,"\n")


#Kompton Kontinuum:

#par, cov=curve_fit(WQ_compton,x_cs[peak[4]-d:peak[4]+d],N_cs[peak[4]-d:peak[4]+d], p0=[12,m*peak[4],1,1] ,sigma=np.sqrt(N_cs[peak[4]-d:peak[4]+d])+0.01)
#par = unc.correlated_values(par, cov)
#h = float(unp.nominal_values(par[0]))
#u = float(unp.nominal_values(par[1]))

plt.figure(constrained_layout=True)
plt.bar(x_cs[peak[2]:peak[5]],N_cs[peak[2]:peak[5]],width=m,label="Messdaten 137Cs")
#plt.plot(m*peak,N_cs[peak],"rx",label="Peak")
plt.xlabel("Energie")
plt.ylabel("Anzahl N")
plt.yscale("log")
#plt.xlim(a*m,b*m)
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
