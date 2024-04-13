import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit 
import uncertainties as unc
import uncertainties.unumpy as unp

x=np.genfromtxt('data1.txt', unpack=True)