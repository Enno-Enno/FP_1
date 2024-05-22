import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
from scipy.signal import find_peaks
import uncertainties as unc
import uncertainties.unumpy as unp
import fktn

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
plt.ylabel(r"Anzahl $N$")
plt.yscale("log")
plt.xlim(a,b)
plt.legend()
plt.savefig("build/plt1_Eu.pdf")
#plt.show()


#Energiezuteilung
Energie=np.array([122,245,344,411,444,779,867,964])
W=np.array([28.58,7.583,26.5,2.234,2.821,12.942,4.245,14.605])/100
x=np.linspace(a,b)

par, cov=curve_fit(fktn.f,peak, Energie)
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
d1 = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
d1_f= float(unp.std_devs(par[1]))

print(m,m_f, d1,d1_f)

plt.figure(constrained_layout=True)
plt.plot(x,fktn.f(x,m,d1),"b-",label="Ausgleichsgerade")
plt.plot(peak,Energie,"rx",label="Peaks")
plt.xlabel("Channel")
plt.ylabel(r"$E \, [\mathrm{KeV}]$")
plt.xlim(a,b)
plt.ylim(0,1000)
plt.legend()
plt.savefig("build/plt2_Fit.pdf")
#plt.show()

print(unc.ufloat(m,m_f)*peak)


#Gaußanpassung
d=25

h=np.zeros(len(peak)); h_f=np.zeros(len(peak))
u=np.zeros(len(peak)); u_f=np.zeros(len(peak))
s=np.zeros(len(peak)); s_f=np.zeros(len(peak))
g=np.zeros(len(peak)); g_f=np.zeros(len(peak))

plt.figure(constrained_layout=True)

for i in range(len(peak)):

    par, cov=curve_fit(fktn.gauß,x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d], p0=[12,peak[i],1,1],sigma=np.sqrt(N_eu[peak[i]-d:peak[i]+d])+0.01)
    par = unc.correlated_values(par, cov)
    h[i] = float(unp.nominal_values(par[0]))
    u[i] = float(unp.nominal_values(par[1]))
    s[i] = float(unp.nominal_values(par[2]))
    g[i] = float(unp.nominal_values(par[3]))
    h_f[i]= float(unp.std_devs(par[0]))
    u_f[i]= float(unp.std_devs(par[1]))
    s_f[i]= float(unp.std_devs(par[2]))
    g_f[i]= float(unp.std_devs(par[3]))


    x=np.linspace(x_eu[peak[i]-d],x_eu[peak[i]+d],1000)
    
    plt.subplot(4,2,i+1)
    #plt.errorbar(x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d],yerr=np.sqrt(N_eu[peak[i]-d:peak[i]+d]),fmt="r")
    plt.bar(x_eu[peak[i]-d:peak[i]+d],N_eu[peak[i]-d:peak[i]+d],width=1,yerr=np.sqrt(N_eu[peak[i]-d:peak[i]+d]),label=f"{i+1}")
    plt.plot(x,fktn.gauß(x,h[i],u[i],s[i],g[i]),"g-")#,label="Gauß-Fit")
    #plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
    plt.xlim(x_eu[peak[i]-15],x_eu[peak[i]+15])
    plt.legend()
    
plt.savefig("build/plt3_Gauß.pdf")
#plt.show()

h=unp.uarray(h,h_f)
u=unp.uarray(u,u_f)
s=unp.uarray(abs(s),s_f)
g=unp.uarray(g,g_f)

print("\n",h, "\n \n", u ,"\n\n" ,s ,"\n\n", g,"\n")

#------------------------------------------------------------------------------------------
#Vollenergie Nachweiß

T=3413
I=np.sqrt(2*np.pi)*h*s
A=unc.ufloat(4130,60)*np.exp(-np.log(2)*(215+23*365)/(13.516*365))
theta=2*np.pi*(1-unp.cos(unp.arctan(22.5/85)))
#theta=2*np.pi*(1-85/np.sqrt(85**2+22.5**2))
print(theta/4/np.pi,A)
print(I)


q=I*4*np.pi/(theta*A*T)/W
print("Q:",q)
Q=unp.nominal_values(q)*100

par, cov=curve_fit(fktn.potenz,m*peak+d1, Q, sigma=unp.std_devs(q),p0=[33,0.8,0])
par = unc.correlated_values(par, cov)
a = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
c = float(unp.nominal_values(par[2]))
a_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))
c_f= float(unp.std_devs(par[2]))

print(a,a_f,b,b_f,c,c_f)

x=np.linspace(0,m*peak[-1]+d1,1000)

plt.figure(constrained_layout=True)
plt.plot(x,fktn.potenz(x,a,b,c),"b-",label="Potenz Fit")
plt.errorbar(m*peak+d1,Q,yerr=unp.std_devs(q)*100,fmt="m.",label="Q Vollenergienachweiswahrscheinlichkeit")
plt.xlabel(r"Energie $E \, [\mathrm{KeV}]$")
plt.ylabel(f"Q [%]")
plt.xlim(0,m*peak[-1]+10)
plt.ylim(0,50)
plt.legend()
plt.savefig("build/plt4_Q.pdf")
#plt.show()

datei = open('build/data.txt','w')
datei.writelines(f"{m} {d1} {a} {b} {c}")