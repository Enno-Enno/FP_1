import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
import fktn

#Aus Plot1.py:
m, d1, a1, b1, c1 = np.genfromtxt("build/data.txt", unpack=True)

#Unbekannt

n_un = np.genfromtxt("000Unbekannt-Spektrum.txt", unpack=True)

a=0;b=8000 #Intervall
x_un=np.zeros(len(n_un[a:b]))
N_un=np.zeros(len(n_un[a:b]))

for index, val in enumerate(n_un[a:b]):
    x_un[index]=(index+a)*m +d1
    N_un[index]=val
    
    
peak,_=find_peaks(N_un,height=1000,distance=50)
peak2,_=find_peaks(N_un[3000:8000],height=100,distance=10)
peak3,_=find_peaks(N_un[5000:8000],height=60,distance=10)
peak2=peak2+3000
peak3=peak3+5000
peak=np.concatenate([peak, peak2, peak3])
print("Peak Energien Unbekannt",peak)
print(N_un[peak])

plt.figure(constrained_layout=True)
plt.bar(x_un,N_un,width=m,label="Messdaten Unbekannt")
plt.plot(peak*m+d1,N_un[peak],"rx",label="Peak")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.ylabel("Anzahl N")
plt.yscale("log")
plt.xlim(a*m,b*m)
plt.legend()
plt.savefig("build/plt9_Un.pdf")
#plt.show()


#Gaußanpassung
d=25

h=np.zeros(len(peak)); h_f=np.zeros(len(peak))
u=np.zeros(len(peak)); u_f=np.zeros(len(peak))
s=np.zeros(len(peak)); s_f=np.zeros(len(peak))
g=np.zeros(len(peak)); g_f=np.zeros(len(peak))

plt.figure(constrained_layout=True)

for i in range(len(peak)):

    par, cov=curve_fit(fktn.gauß,x_un[peak[i]-d:peak[i]+d],N_un[peak[i]-d:peak[i]+d], p0=[12,m*peak[i]+d1,1,1],sigma=np.sqrt(N_un[peak[i]-d:peak[i]+d])+0.01)
    par = unc.correlated_values(par, cov)
    h[i] = float(unp.nominal_values(par[0]))
    u[i] = float(unp.nominal_values(par[1]))
    s[i] = float(unp.nominal_values(par[2]))
    g[i] = float(unp.nominal_values(par[3]))
    h_f[i]= float(unp.std_devs(par[0]))
    u_f[i]= float(unp.std_devs(par[1]))
    s_f[i]= float(unp.std_devs(par[2]))
    g_f[i]= float(unp.std_devs(par[3]))

    x=np.linspace(x_un[peak[i]-d],x_un[peak[i]+d],1000)
    
    plt.subplot(7,2,i+1)
    plt.bar(x_un[peak[i]-d:peak[i]+d],N_un[peak[i]-d:peak[i]+d],width=m,yerr=np.sqrt(N_un[peak[i]-d:peak[i]+d]),label=f"Messdaten Peak {i+1}")
    plt.plot(x,fktn.gauß(x,h[i],u[i],s[i],g[i]),"g-",label="Gauß-Fit")
    plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
    plt.xlim(x_un[peak[i]-d],x_un[peak[i]+d])
    plt.legend()
    
plt.show()

h=unp.uarray(h,h_f)
u=unp.uarray(u,u_f)
s=unp.uarray(abs(s),s_f)
g=unp.uarray(g,g_f)

I=np.sqrt(2*np.pi)*h*s
print(I)
print(u)