import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

d, B=np.genfromtxt('data1.txt', unpack=True)

th1_n1,th1_n1_m,th2_n1,th2_n1_m, th1_r,th1_r_m,th2_r,th2_r_m, th1_n2,th1_n2_m,th2_n2,th2_n2_m=np.genfromtxt('data1.txt', unpack=True)



plt.plot(d,B, "rx",label="Messdaten B-Feld")
plt.legend()
plt.show()