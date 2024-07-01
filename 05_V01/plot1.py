import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

def f(x, m, d):
    return m*x+d 



n_eu=np.genfromtxt('data1.txt', unpack=True)

x=np.zeros(len(n_eu))
N=np.zeros(len(n_eu))

for index, val in enumerate(n_eu):
    x[index]=index
    N[index]=val

plt.plot(x,np.log(N),"kx")
plt.show()

par, cov=curve_fit(f,x, np.log(N))
par = unc.correlated_values(par, cov)
m = float(unp.nominal_values(par[0]))
d1 = float(unp.nominal_values(par[1]))
m_f= float(unp.std_devs(par[0]))
d1_f= float(unp.std_devs(par[1]))
x_p=np.linspace(0,100)

plt.plot(x,np.log(N),"kx")
plt.plot(x_p,f(x_p,m,d1),"g-")
plt.show()

print(m,m_f, d1,d1_f)