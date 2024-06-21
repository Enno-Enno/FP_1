import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import scipy.constants as con

def f(x, m):
    return m*x#+b 

d, B=np.genfromtxt('data1.txt', unpack=True)
th1_n1,th1_n1_m,th2_n1,th2_n1_m, th1_r,th1_r_m,th2_r,th2_r_m, th1_n2,th1_n2_m,th2_n2,th2_n2_m=np.genfromtxt('data2.txt', unpack=True)

th1_n1=th1_n1+th1_n1_m/60
th2_n1=th2_n1+th2_n1_m/60
th1_r=th1_r+th1_r_m/60
th2_r=th2_r+th2_r_m/60
th1_n2=th1_n2+th1_n2_m/60
th2_n2=th2_n2+th2_n2_m/60

L_n1=1.296*10**(-3)
L_r=5.11*10**(-3)
L_n2=1.36*10**(-3)

lamda=np.array([1.06,1.29,1.45,1.72,1.96,2.156,2.34,2.51,2.65])#*10**(-6)
d_x=np.linspace(68,132)
x=np.linspace(1,8)#*10**(-12)

delta_n1=abs(th1_n1-th2_n1)/(2*L_n1)/180*np.pi
delta_r=abs(th1_r-th2_r)/(2*L_r)/180*np.pi
delta_n2=abs(th1_n2-th2_n2)/(2*L_n2)/180*np.pi

print(delta_n1)
print(delta_r)
print(delta_n2)

plt.figure(constrained_layout=True)
plt.plot(d,B, "rx",label="Messdaten B-Feld")
plt.plot(d_x,0.000000001*x+428,"b--",label="Maximales Feld")
plt.xlabel(r"Abstand $d \, [\mathrm{mm}]$")
plt.ylabel(r"$B \, [\mathrm{mT}]$")
plt.xlim(d_x[0],d_x[-1])
plt.legend()
plt.savefig("build/B_Feld.pdf")
#plt.show()

plt.figure(constrained_layout=True)
plt.plot(lamda**2,delta_n1, "rx",label="Messdaten Probe 1")
plt.plot(lamda**2,delta_r, "bx",label="Messdaten Probe 2")
plt.plot(lamda**2,delta_n2, "gx",label="Messdaten Probe 3")
plt.ylabel(r"$\theta  [\mathrm{rad/m}]$")
plt.xlabel(r"$\lambda^2 [\mathrm{\mu m}^2]$")
plt.legend()
plt.savefig("build/Messdaten.pdf")
#plt.show()

par, cov=curve_fit(f,lamda**2, abs(delta_n2-delta_r))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
m_f= float(unp.std_devs(par[0]))
print(m,m_f)

plt.figure(constrained_layout=True)
plt.plot(lamda**2,abs(delta_n2-delta_r), "gx",label=r"$\Theta_{2,3}$")
plt.plot(x,f(x,m),"b-", label="Ausgleichsgerade")
plt.ylabel(r"$[\mathrm{rad/m}]$")
plt.xlabel(r"$\lambda^2 [\mathrm{\mu m}^2]$")
plt.legend()
plt.savefig("build/plot1.pdf")
#plt.show()

m=unp.uarray(m,m_f)
m1=unp.sqrt((con.e**3*1.28*10**(24)*0.428)/(8*3.857*con.epsilon_0*m*10**(12)*con.c**3*np.pi**2))

par, cov=curve_fit(f,lamda**2, abs(delta_n1-delta_r))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
m_f= float(unp.std_devs(par[0]))
print(m,m_f)

plt.figure(constrained_layout=True)
plt.plot(lamda**2,abs(delta_n1-delta_r), "gx",label=r"$\Theta_{2,1}$")
plt.plot(x,f(x,m),"b-", label="Ausgleichsgerade")
plt.ylabel(r"$[\mathrm{rad/m}]$")
plt.xlabel(r"$\lambda^2 [\mathrm{\mu m}^2]$")
plt.legend()
plt.savefig("build/plot2.pdf")
#plt.show()

m=unp.uarray(m,m_f)
m2=unp.sqrt((con.e**3*2.8*10**(24)*0.428)/(8*3.857*con.epsilon_0*m*10**(12)*con.c**3*np.pi**2))

print(m1/con.m_e)
print(m2/con.m_e)
print((m1+m2)/2/con.m_e)