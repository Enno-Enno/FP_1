import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

def f(x, m, b):
    return m*x+b 

d, B=np.genfromtxt('data1.txt', unpack=True)
th1_n1,th1_n1_m,th2_n1,th2_n1_m, th1_r,th1_r_m,th2_r,th2_r_m, th1_n2,th1_n2_m,th2_n2,th2_n2_m=np.genfromtxt('data2.txt', unpack=True)

th1_n1=th1_n1+th1_n1_m/60
th2_n1=th2_n1+th2_n1_m/60
th1_r=th1_r+th1_r_m/60
th2_r=th2_r+th2_r_m/60
th1_n2=th1_n2+th1_n2_m/60
th2_n2=th2_n2+th2_n2_m/60

L_n1=1.296
L_r=5.11
L_n2=1.36

lamda=np.array([1.06,1.29,1.45,1.72,1.96,2.156,2.34,2.51,2.65])
x=np.linspace(1,8)

delta_n1=abs(th1_n1-th2_n1)/(2*L_n1)/180*np.pi
delta_r=abs(th1_r-th2_r)/(2*L_r)/180*np.pi
delta_n2=abs(th1_n2-th2_n2)/(2*L_n2)/180*np.pi

print(delta_n1)

plt.figure(constrained_layout=True)
plt.plot(d,B, "rx",label="Messdaten B-Feld")
plt.legend()
plt.show()

plt.figure(constrained_layout=True)
plt.plot(lamda**2,delta_n1, "rx",label="Messdaten Probe 1")
plt.ylabel(r"$\theta  [\mathrm{°/mm}]$")
plt.xlabel(r"$\lambda [\mathrm{\mu m}^2]$")
plt.legend()
plt.show()

plt.figure(constrained_layout=True)
plt.plot(lamda**2,delta_r, "rx",label="Messdaten Probe 2")
plt.ylabel(r"$\theta  [\mathrm{°/mm}]$")
plt.xlabel(r"$\lambda [\mathrm{\mu m}^2]$")
plt.legend()
plt.show()

plt.figure(constrained_layout=True)
plt.plot(lamda**2,delta_n2, "rx",label="Messdaten Probe 3")
plt.ylabel(r"$\theta  [\mathrm{°/mm}]$")
plt.xlabel(r"$\lambda [\mathrm{\mu m}^2]$")
plt.legend()
plt.show()

par, cov=curve_fit(f,lamda**2, abs(delta_n2-delta_r))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

print(m,m_f)
print(b,b_f)

plt.figure(constrained_layout=True)
plt.plot(lamda**2,abs(delta_n2-delta_r), "rx",label="Weiß noch nicht")
plt.plot(x,f(x,m,b),"b-", label="Ausgleichsgerade")
plt.ylabel(r"$\theta  [\mathrm{°/mm}]$")
plt.xlabel(r"$\lambda [\mathrm{\mu m}^2]$")
plt.legend()
plt.show()

par, cov=curve_fit(f,lamda**2, abs(delta_n1-delta_r))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

print(m,m_f)
print(b,b_f)

plt.figure(constrained_layout=True)
plt.plot(lamda**2,abs(delta_n1-delta_r), "rx",label="Weiß noch nicht")
plt.plot(x,f(x,m,b),"b-", label="Ausgleichsgerade")
plt.ylabel(r"$\theta  [\mathrm{°/mm}]$")
plt.xlabel(r"$\lambda [\mathrm{\mu m}^2]$")
plt.legend()
plt.show()