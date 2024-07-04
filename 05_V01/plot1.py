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


N_kal=np.genfromtxt('data_kal.txt', unpack=True)
N_myon=np.genfromtxt('data1.txt', unpack=True)
T,N=np.genfromtxt('data2.txt', unpack=True)

#Plato Bestimmung 

par, cov=curve_fit(f,T[3:9], N[3:9])
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

print(m,b)
x=np.linspace(4,16)

plt.plot(x,f(x,m,b),"r-",label="Ausgleichsgerade")
plt.plot(T[3:9],N[3:9],"bx",label="Messwerte Plato")
plt.plot(T[9:-1],N[9:-1],"gx",label="Messwerte Randbereich")
plt.plot(T[0:3],N[0:3],"gx")
plt.ylabel(r"$N$")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.show()


#Kalibrierung

peak,_=find_peaks(N_kal,distance=10,height=10)
t_kal=np.arange(1, len(peak)+1)*5-1
t_kal=t_kal/10
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
plt.ylabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.show()


#Myon Lebenszeitbestimmung (hoffentlich)

t_myon=np.arange(1, len(N_myon)+1)
t_myon=f(t_myon,m,b)

A=int(np.floor((1.1-b)/m)+1)
B=int(np.floor((10-b)/m)+1)
print(A,B)

plt.bar(t_myon[A:B],N_myon[A:B],width=m,label="Messwerte")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.show()

par, cov=curve_fit(exp,t_myon[A:B],N_myon[A:B],p0=[1/2.2,125,2])#,sigma=np.sqrt(N_myon[A:B]))
par = unc.correlated_values(par, cov)
lamda = float(unp.nominal_values(par[0]))
N = float(unp.nominal_values(par[1]))
U = float(unp.nominal_values(par[2]))
lamda_f= float(unp.std_devs(par[0]))
N_f= float(unp.std_devs(par[1]))
U_f= float(unp.std_devs(par[2]))

par, cov=curve_fit(f,t_myon[A:B], np.log(N_myon[A:B]),sigma=np.sqrt(N_myon[A:B]))
par = unc.correlated_values(par, cov)
m2 = float(unp.nominal_values(par[0]))
b2 = float(unp.nominal_values(par[1]))
m2_f= float(unp.std_devs(par[0]))
b2_f= float(unp.std_devs(par[1]))

t=np.linspace(A*m+b,B*m+b)

plt.bar(t_myon[A:B],np.log(N_myon[A:B]),width=m,label="Messwerte")
plt.plot(t,f(t,m2,b2),"g-",label="Ausgleichsgerade")
plt.plot(t,np.log(exp(t,lamda,N,U)),"r-")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.show()

lamda=unp.uarray(lamda,lamda_f)
m2=unp.uarray(m2,m2_f)
print(1/lamda)
print(1/m2)
print(lamda,N,U)