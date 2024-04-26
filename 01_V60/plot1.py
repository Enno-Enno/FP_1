import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit 
import uncertainties as unc
import uncertainties.unumpy as unp

x=np.genfromtxt('data1.txt', unpack=True)

Schwellwert_unten = np.array([34.8,34.6])
Schwellwert_oben = np.array([36,34.9, 34.8])