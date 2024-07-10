import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks

def g(x, c):
    return 0*x+c

def f(x, m, b):
    return m*x+b 

def exp(t,lamda,N,U):
    return N*np.exp(-lamda*t) +U

def poisson(t,lamda,N,U):
    return N*lamda*t*np.exp(-lamda*t) +U


N_kal=np.genfromtxt('data_kal.txt', unpack=True)
N_myon=np.genfromtxt('data1.txt', unpack=True)
T,N=np.genfromtxt('data2.txt', unpack=True)

#Plato Bestimmung 

par, cov=curve_fit(g,T[3:9], N[3:9],sigma=np.sqrt(N[3:9]))
par = unc.correlated_values(par, cov)
c = float(unp.nominal_values(par[0]))
#b = float(unp.nominal_values(par[1]))
c_f= float(unp.std_devs(par[0]))
#b_f= float(unp.std_devs(par[1]))

print(c)
x=np.linspace(4,14)

plt.figure(constrained_layout=True)
plt.plot(x,g(x,c),"b-",label="Plateau")
plt.errorbar(T[3:9],N[3:9],yerr=np.sqrt(N[3:9]),fmt="gx",label="Messwerte Plateau")
plt.errorbar(T[9:-1],N[9:-1],yerr=np.sqrt(N[9:-1]),fmt="rx",label="Messwerte Randbereich")
plt.errorbar(T[0:3],N[0:3],yerr=np.sqrt(N[0:3]),fmt="rx")
plt.ylabel(r"$N$")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.savefig("build/plot1.pdf")
#plt.show()


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

plt.figure(constrained_layout=True)
plt.plot(x,f(x,m,b),"b-",label="Ausgleichsgerade")
plt.plot(peak,t_kal,"gx",label="Messwerte Kalibrierung")
plt.xlabel(r"Channel")
plt.ylabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.savefig("build/plot2.pdf")
#plt.show()


#Myon Lebenszeitbestimmung (hoffentlich)

t_myon=np.arange(1, len(N_myon)+1)
t_myon=f(t_myon,m,b)

A=int(np.floor((1.1-b)/m)+1)-2
B=int(np.floor((10-b)/m)+1)
C=17
print(A,B)

plt.figure(constrained_layout=True)
plt.bar(t_myon[A:B],np.log(N_myon[A:B]),color="c",width=m,label="Messwerte NACH der Kante")
plt.bar(t_myon[C:A],np.log(N_myon[C:A]),color="b",width=m,label="Messwerte VOR der Kante")
plt.bar(t_myon[0:C],np.log(N_myon[0:C]),color="r",width=m,label="ungenutzte Werte")
plt.bar(t_myon[B:-1],np.log(N_myon[B:-1]),color="r",width=m)
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"Anzahl $\log(N)$")
plt.xlim(0,t_myon[-1])
plt.legend()
plt.savefig("build/plot3.pdf")
#plt.show()

par, cov=curve_fit(exp,t_myon[A:B],N_myon[A:B],p0=[1/2.2,125,2],sigma=np.sqrt(N_myon[A:B]))
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

t=np.linspace(A*m+b,B*m+b,1000)

plt.figure(constrained_layout=True)
plt.bar(t_myon[A:B],np.log(N_myon[A:B]),color="c",width=m,label="Messwerte")
plt.plot(t,f(t,m2,b2),"b-",label="Ausgleichsgerade")
plt.plot(t,np.log(exp(t,lamda,N,U)),"r-",label="Exponitialfit")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.xlim(A*m+b,B*m+b)
plt.legend()
plt.savefig("build/plot4.pdf")
#plt.show()

lamda=unp.uarray(lamda,lamda_f)
m2=unp.uarray(m2,m2_f)
print("Tau=",1/lamda)
print("Tau=",1/m2)
print(lamda,N,U)
print(m2,b2)

#Poisson Verteilung #Kann man eigentlich ignorieren
t=np.linspace(0,B*m+b,1000)

par, cov=curve_fit(poisson,t_myon[3:B],np.log(N_myon[3:B]),p0=[2.2,50000,0])#,sigma=np.sqrt(np.log(N_myon[3:B])))
par = unc.correlated_values(par, cov)
lamda = float(unp.nominal_values(par[0]))
N = float(unp.nominal_values(par[1]))
U = float(unp.nominal_values(par[2]))
lamda_f= float(unp.std_devs(par[0]))
N_f= float(unp.std_devs(par[1]))
U_f= float(unp.std_devs(par[2]))

plt.figure(constrained_layout=True)
plt.bar(t_myon[3:B],np.log(N_myon[3:B]),color="c",width=m,label="Messwerte")
plt.plot(t,poisson(t,lamda,N,U),"r-",label="Poissonfit")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"Anzahl $\log(N)$")
plt.xlim(0,B*m+b)
plt.legend()
#plt.show()

print(lamda,N,U)
print("Tau=",1/lamda)

#Anpassung der Kannte
N_myon[C:A]=N_myon[C:A]-(N_myon[A-3]-N_myon[A])

par, cov=curve_fit(exp,t_myon[C:B],N_myon[C:B],p0=[1/2.2,125,2],sigma=np.sqrt(N_myon[C:B]))
par = unc.correlated_values(par, cov)
lamda = float(unp.nominal_values(par[0]))
N = float(unp.nominal_values(par[1]))
U = float(unp.nominal_values(par[2]))
lamda_f= float(unp.std_devs(par[0]))
N_f= float(unp.std_devs(par[1]))
U_f= float(unp.std_devs(par[2]))

t=np.linspace(C*m+b,B*m+b,1000)

plt.figure(constrained_layout=True)
plt.bar(t_myon[A:B],np.log(N_myon[A:B]),color="c",width=m,label="Messwerte NACH der Kante")
plt.bar(t_myon[C:A],np.log(N_myon[C:A]),color="b",width=m,label="Korregierte Messwerte VOR der Kante")
plt.plot(t,np.log(exp(t,lamda,N,U)),"r-",label="Exponitialfit")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"Anzahl $\log(N)$")
plt.xlim(C*m+b,B*m+b)
plt.legend()
plt.savefig("build/plot5.pdf")
#plt.show()

print("Tau=",1/lamda)
print(lamda,N,U)

