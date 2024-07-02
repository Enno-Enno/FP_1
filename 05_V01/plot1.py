import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks

def f(x, m, b):
    return m*x+b 

def exp(t,lamda,N,U):
    return N*lamda*t*np.exp(-lamda*t) +U


N_kal=np.genfromtxt('data2.txt', unpack=True)
N_myon=np.genfromtxt('data1.txt', unpack=True)

peak,_=find_peaks(N_kal,distance=10,height=10)
t_kal=np.arange(1, len(peak)+1)*5-1
print(peak)
print(t_kal)

par, cov=curve_fit(f,peak, t_kal)
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

print(m,b)

x=np.linspace(0,450)

plt.plot(x,f(x,m,b),"m-",label="Ausgleichsgerade")
plt.plot(peak,t_kal,"bx",label="Messwerte Kalibrierung")
plt.xlabel(r"Channel")
plt.ylabel(r"$t [\mathrm{\mu s}]$")
plt.legend()
plt.show()

t_myon=np.arange(1, len(N_myon)+1)
t_myon=f(t_myon,m,b)

plt.plot(t_myon,np.log(N_myon),"bx",label="Messwerte")
plt.xlabel(r"$t [\mathrm{\mu s}]$")
plt.legend()
plt.show()


A=17;B=int(np.floor((10-b)/m)+1)
print(B)

#par, cov=curve_fit(exp,t_myon[A:B],N_myon[A:B])#,sigma=np.sqrt(N_myon[A:B]),p0=[2,125281093,0])
#par = unc.correlated_values(par, cov)
#lamda = float(unp.nominal_values(par[0]))
#N = float(unp.nominal_values(par[1]))
#U = float(unp.nominal_values(par[2]))
#lamda_f= float(unp.std_devs(par[0]))
#N_f= float(unp.std_devs(par[1]))
#U_f= float(unp.std_devs(par[2]))

par, cov=curve_fit(f,t_myon[A:B], np.log(N_myon[A:B]))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

t=np.linspace(2,10)

plt.plot(t,f(t,m,b),"m-",label="Fit")
plt.plot(t_myon[A:B],np.log(N_myon[A:B]),"bx",label="Messwerte")
plt.xlabel(r"$t [\mathrm{\mu s}]$")
#plt.yscale("log")
plt.xlim(0,60)
plt.legend()
plt.show()

print(np.log(2)/m)