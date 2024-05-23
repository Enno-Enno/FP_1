import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
from scipy.integrate import quad
import fktn ## Enthält die Hilfreichen Funktionen fürs Compton Fitten :)


#Aus Plot1.py:
m,_, d1,_, a1, b1, c1 = np.genfromtxt("build/data.txt", unpack=True)

#Funktionen:
## In anderer Datei!

WQ_compton = lambda E_gamma,theta: fktn.WQ_compton


n_cs = np.genfromtxt("137Cs-Spektrum.txt", unpack=True)

a=0;b=4000 #Intervall
x_cs=np.zeros(len(n_cs[a:b]))
N_cs=np.zeros(len(n_cs[a:b]))

for index, val in enumerate(n_cs[a:b]):
    x_cs[index]=(index+a)*m +d1
    N_cs[index]=val
    
    
peak,_=find_peaks(N_cs,height=70,distance=50)
peak2,_=find_peaks(N_cs[2200:2500],height=42,distance=50)
peak2=peak2+2200
peak=np.concatenate([peak, peak2])
print("Peak Energien Cs",peak)
print(N_cs[peak])


plt.figure(constrained_layout=True)
plt.bar(x_cs,N_cs,width=m,label="Messdaten 137Cs")
plt.plot(m*peak+d1,N_cs[peak],"rx",label="Peak")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
plt.savefig("build/plt5_Cs.pdf")
#plt.show()

#Zentelbreite und Halbwertsbreite Inhalt Gauß
d=15

par, cov=curve_fit(fktn.gauß,x_cs[peak[2]-d:peak[2]+d],N_cs[peak[2]-d:peak[2]+d], p0=[12,m*peak[2]+d1,1,1] ,sigma=np.sqrt(N_cs[peak[2]-d:peak[2]+d])+0.01)
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

delta_xZ=s*np.sqrt(2*np.log(10*h/(h-9*g))) #Zentelbreite
xZ=np.linspace(u-delta_xZ,u+delta_xZ,1000)


x=np.linspace(x_cs[peak[2]-d],x_cs[peak[2]+d],1000)

plt.figure(constrained_layout=True)
plt.bar(x_cs[peak[2]-d:peak[2]+d],N_cs[peak[2]-d:peak[2]+d],width=m,yerr=np.sqrt(N_cs[peak[2]-d:peak[2]+d]),label=f"Messdaten Photopeak ")
plt.plot(x,fktn.gauß(x,h,u,s,g),"g-",label="Gauß-Fit")
plt.plot(xH,fktn.f(xH,0.000001,0.5*h),"y--", label="Halbwertsbreite")
plt.plot(xZ,fktn.f(xZ,0.000001,0.1*h),"r--", label="Zehntelwertsbreite")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.xlim(u-3,u+3)
plt.legend()
plt.savefig("build/plt6_Ph_peak.pdf")
#plt.show()

h=unp.uarray(h,h_f)
u=unp.uarray(u,u_f)
s=unp.uarray(abs(s),s_f)
g=unp.uarray(g,g_f)

print("Zentelbreite=",2*s*unp.sqrt(2*unp.log(10*h/(h-9*g))))
print("Halbwertsbreite=",2*s*unp.sqrt(2*unp.log(2*h/(h-g))))
print("Inhalt=",np.sqrt(2*np.pi)*h*s/m)

print("\n",h, "\n ", u ,"\n" ,s ,"\n", g,"\n")


#Kompton Kontinuum:

E_Cs = 661.454
Comtonkante = fktn.E_kante(E_Cs)
Rückstrahlpeak = fktn.E_rück(E_Cs)
d=50

print(Comtonkante,Rückstrahlpeak)

par, cov=curve_fit(fktn.WQ_Energie,x_cs[peak[1]+d:peak[-1]+d],N_cs[peak[1]+d:peak[-1]+d],p0=[E_Cs,3,10],sigma=np.sqrt(N_cs[peak[1]+d:peak[-1]+d]))
par = unc.correlated_values(par, cov)
E = float(unp.nominal_values(par[0]))
k = float(unp.nominal_values(par[1]))
unterG = float(unp.nominal_values(par[2]))
E_f = float(unp.std_devs(par[0]))
k_f = float(unp.std_devs(par[1]))
unterG_f = float(unp.std_devs(par[2]))

print(E,k,unterG)
x=np.linspace(Rückstrahlpeak,Comtonkante,1000)

plt.figure(constrained_layout=True)
plt.bar(x_cs[peak[1]-d:peak[-1]+150],N_cs[peak[1]-d:peak[-1]+150],width=m,label="Messdaten 137Cs")
plt.plot(x,fktn.WQ_Energie(x,E,k,unterG),"m", label= "Compton fit ")
plt.axvline(Comtonkante,color="m",ls="--",label="Comtonkante")
plt.axvline(Rückstrahlpeak,color="m",ls="-.",label="Rückstrahlpeak")
plt.plot(m*peak+d1,N_cs[peak],"rx",label="Peak")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(x_cs[peak[1]-d],x_cs[peak[-1]+150])
plt.ylim(1,200)
plt.legend()
plt.savefig("build/plt7_Compton.pdf")
#plt.show()

Z_com,_ = quad(fktn.WQ_Energie,x_cs[peak[1]+d], x_cs[peak[-1]+d], args=(E,k,unterG))
print(Z_com/m, np.sqrt(Z_com/m))
print(Z_com/sum(N_cs))

#Absorbtionswahrscheinlichkeit
l=3.9
u_c=0.437
u_ph=0.369363

p_c=1-np.exp(-u_c*l)
p_ph=1-np.exp(-u_ph*l)
print(p_c,p_ph)

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#Aktivitätsbestimmung

n_ba = np.genfromtxt("133Ba-Spektrum.txt", unpack=True)

W=np.array([34.06,7.164,18.33,62.05,8.94])/100

a=0;b=3000 #Intervall
x_ba=np.zeros(len(n_ba[a:b]))
N_ba=np.zeros(len(n_ba[a:b]))

for index, val in enumerate(n_ba[a:b]):
    x_ba[index]=(index+a)*m+d1
    N_ba[index]=val

    
peak_ba,_=find_peaks(N_ba,height=100)
print("Peak Energien Ba",m*peak_ba+d1)
print(N_ba[peak_ba])

plt.figure(constrained_layout=True)
plt.bar(x_ba,N_ba,width=m,label="Messdaten 133Ba")
plt.plot(m*peak_ba+d1,N_ba[peak_ba],"rx",label="Peak")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
plt.savefig("build/plt8_Ba.pdf")
#plt.show()

#Gaußanpassung
d=25

x_ba=(x_ba-d1)/m

h=np.zeros(len(peak_ba)); h_f=np.zeros(len(peak_ba))
u=np.zeros(len(peak_ba)); u_f=np.zeros(len(peak_ba))
s=np.zeros(len(peak_ba)); s_f=np.zeros(len(peak_ba))
g=np.zeros(len(peak_ba)); g_f=np.zeros(len(peak_ba))

plt.figure(constrained_layout=True)

for i in range(len(peak_ba)):

    par, cov=curve_fit(fktn.gauß,x_ba[peak_ba[i]-d:peak_ba[i]+d],N_ba[peak_ba[i]-d:peak_ba[i]+d], p0=[12,peak_ba[i],1,1],sigma=np.sqrt(N_ba[peak_ba[i]-d:peak_ba[i]+d])+0.01)
    par = unc.correlated_values(par, cov)
    h[i] = float(unp.nominal_values(par[0]))
    u[i] = float(unp.nominal_values(par[1]))
    s[i] = float(unp.nominal_values(par[2]))
    g[i] = float(unp.nominal_values(par[3]))
    h_f[i]= float(unp.std_devs(par[0]))
    u_f[i]= float(unp.std_devs(par[1]))
    s_f[i]= float(unp.std_devs(par[2]))
    g_f[i]= float(unp.std_devs(par[3]))

    x=np.linspace(x_ba[peak_ba[i]-d],x_ba[peak_ba[i]+d],1000)
    
    plt.subplot(3,2,i+1)
    plt.bar(x_ba[peak_ba[i]-d:peak_ba[i]+d],N_ba[peak_ba[i]-d:peak_ba[i]+d],width=1,yerr=np.sqrt(N_ba[peak_ba[i]-d:peak_ba[i]+d]),label=f"Messdaten Peak {i+1}")
    plt.plot(x,fktn.gauß(x,h[i],u[i],s[i],g[i]),"g-",label="Gauß-Fit")
    plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
    plt.xlim(x_ba[peak_ba[i]-d],x_ba[peak_ba[i]+d])
    plt.legend()
    
#plt.savefig("build/plt3_Gauß.pdf")
#plt.show()

h=unp.uarray(h,h_f)
u=unp.uarray(u,u_f)
s=unp.uarray(abs(s),s_f)
g=unp.uarray(g,g_f)

I=np.sqrt(2*np.pi)*h*s
theta=2*np.pi*(1-unp.cos(unp.arctan(22.5/85)))
T=3816
Q=fktn.potenz(peak_ba*m+d1,a1,b1,c1)/100
A=I*4*np.pi/(theta*Q*W*T)

print("I=",I)
print(Q)
print(A)
print(np.mean(A))

#------------------------------------------------------------------------------------------------