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
c_f= float(unp.std_devs(par[0]))

print(c)
x=np.linspace(3,15)

plt.figure(constrained_layout=True)
plt.plot(x,g(x,c),"b-",label="Plateau")
plt.errorbar(T[3:9],N[3:9],yerr=np.sqrt(N[3:9]),fmt="gx",label="Messwerte Plateau")
plt.errorbar(T[9:],N[9:],yerr=np.sqrt(N[9:]),fmt="rx",label="Messwerte Randbereich")
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

plt.figure(constrained_layout=True)
plt.bar(np.arange(1, len(N_kal)+1),N_kal,color="b",label="Messwerte")
plt.xlabel(r"Channel")
plt.ylabel(r"$N$")
plt.yscale("log")
plt.xlim(0,460)
plt.legend()
plt.savefig("build/plot2.pdf")

par, cov=curve_fit(f,peak, t_kal)
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
b = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
b_f= float(unp.std_devs(par[1]))

print(m,m_f,b,b_f)

x=np.linspace(0,450)

plt.figure(constrained_layout=True)
plt.plot(x,f(x,m,b),"b-",label="Ausgleichsgerade")
plt.plot(peak,t_kal,"gx",label="Messwerte Kalibrierung")
plt.xlabel(r"Channel")
plt.ylabel(r"$t \, [\mathrm{\mu s}]$")
plt.legend()
plt.savefig("build/plot3.pdf")
#plt.show()


#Myon Lebenszeitbestimmung (hoffentlich)

t_myon=np.arange(1, len(N_myon)+1)
t_myon=f(t_myon,m,b)

A=int(np.floor((1.1-b)/m)+1)-2
B=int(np.floor((10-b)/m)+1)-1
C=17
print(A,B,C)
print(m*np.array([A,B,C])+b)
print(sum(N_myon))

plt.figure(constrained_layout=True)
plt.bar(t_myon[A:B],N_myon[A:B],color="c",width=m,label="Messwerte NACH der Kante")
plt.bar(t_myon[C:A],N_myon[C:A],color="b",width=m,label="Messwerte VOR der Kante")
plt.bar(t_myon[0:C],N_myon[0:C],color="r",width=m,label="ungenutzte Werte")
plt.bar(t_myon[B:-1],N_myon[B:-1],color="r",width=m)
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"Anzahl $N$")
plt.xlim(0,t_myon[-1])
plt.yscale("log")
plt.legend()
plt.savefig("build/plot4.pdf")
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
plt.bar(t_myon[A:B],N_myon[A:B],color="c",width=m,label="Messwerte")#,yerr=np.sqrt(N_myon[A:B])
plt.plot(t,exp(t,lamda,N,U),"k-",label="Exponitialfit")
plt.plot(t,np.exp(f(t,m2,b2)),"m-",label="Ausgleichsgerade")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"$N$")
plt.xlim(A*m+b,B*m+b)
plt.yscale("log")
plt.legend()
plt.savefig("build/plot5.pdf")
#plt.show()

lamda=unp.uarray(lamda,lamda_f)
N=unp.uarray(N,N_f)
U=unp.uarray(U,U_f)
m2=unp.uarray(m2,m2_f)
b2=unp.uarray(b2,b2_f)
print("Tau=",1/lamda)
print("Tau=",1/m2)
print(lamda,N,U)
print(m2,b2,unp.exp(b2))

#Poisson Verteilung #Kann man eigentlich ignorieren
#t=np.linspace(0,B*m+b,1000)
#
#par, cov=curve_fit(poisson,t_myon[3:B],np.log(N_myon[3:B]),p0=[2.2,50000,0])#,sigma=np.sqrt(np.log(N_myon[3:B])))
#par = unc.correlated_values(par, cov)
#lamda = float(unp.nominal_values(par[0]))
#N = float(unp.nominal_values(par[1]))
#U = float(unp.nominal_values(par[2]))
#lamda_f= float(unp.std_devs(par[0]))
#N_f= float(unp.std_devs(par[1]))
#U_f= float(unp.std_devs(par[2]))
#
#plt.figure(constrained_layout=True)
#plt.bar(t_myon[3:B],np.log(N_myon[3:B]),color="c",width=m,label="Messwerte")
#plt.plot(t,poisson(t,lamda,N,U),"r-",label="Poissonfit")
#plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
#plt.ylabel(r"$\log(N)$")
#plt.xlim(0,B*m+b)
#plt.legend()
##plt.show()
#
#print(lamda,N,U)
#print("Tau=",1/lamda)

#Anpassung der Kannte
#B=383
N_myon[C:A]=N_myon[C:A]-(N_myon[A-3]-N_myon[A])-100

par, cov=curve_fit(exp,t_myon[C:B],N_myon[C:B],p0=[1/2.2,20,2],sigma=np.sqrt(N_myon[C:B]))
par = unc.correlated_values(par, cov)
lamda2 = float(unp.nominal_values(par[0]))
N2 = float(unp.nominal_values(par[1]))
U2 = float(unp.nominal_values(par[2]))
lamda2_f= float(unp.std_devs(par[0]))
N2_f= float(unp.std_devs(par[1]))
U2_f= float(unp.std_devs(par[2]))

par, cov=curve_fit(f,t_myon[C:B], np.log(N_myon[C:B]),sigma=np.sqrt(N_myon[C:B]))
par = unc.correlated_values(par, cov)
m3 = float(unp.nominal_values(par[0]))
b3 = float(unp.nominal_values(par[1]))
m3_f= float(unp.std_devs(par[0]))
b3_f= float(unp.std_devs(par[1]))

t=np.linspace(C*m+b,B*m+b,1000)

plt.figure(constrained_layout=True)
plt.bar(t_myon[A:B],N_myon[A:B],color="c",width=m,label="Messwerte NACH der Kante")
plt.bar(t_myon[C:A],N_myon[C:A],color="b",width=m,label="Korregierte Messwerte VOR der Kante")
plt.plot(t,exp(t,lamda2,N2,U2),"k-",label="Exponitialfit")
plt.plot(t,np.exp(f(t,m3,b3)),"m-",label="Ausgleichsgerade")
plt.xlabel(r"$t \, [\mathrm{\mu s}]$")
plt.ylabel(r"$N$")
plt.xlim(C*m+b,B*m+b)
plt.yscale("log")
plt.legend()
plt.savefig("build/plot6.pdf")
#plt.show()

lamda2=unp.uarray(lamda2,lamda2_f)
N2=unp.uarray(N2,N2_f)
U2=unp.uarray(U2,U2_f)
m3=unp.uarray(m3,m3_f)
b3=unp.uarray(b3,b3_f)
print("Tau=",1/lamda2)
print("Tau=",1/m3)
print(lamda2,N2,U2)
print(m3,b3,unp.exp(b3))

print("schummeln:",(abs(1/m2+1/m3)+1/lamda)/3)