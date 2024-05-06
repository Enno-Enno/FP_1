import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import uncertainties as unc
import uncertainties.unumpy as unp

n = np.genfromtxt("EU-Spektrum.txt", unpack=True)

def WQ_compton(alpha,theta):
    r_0, _, delta_r = const["classical electron radius"]
    Z = 40
    cos = lambda x: np.cos(x)
    WQ = (
        Z
        * r_0**2
        * (1 / (alpha * (1 - cos(theta)) + 1)) ** 2
        * (cos(theta) ** 2 + 1) / 2
        * (
            (alpha**2 * (1 - cos(theta)) ** 2)
            / (((alpha * (-cos(theta) + 1) + 1) * (cos(theta) ** 2 + 1)))
            + 1
        )
    )  # - Derivative(sigma, Omega) == 0
    return WQ
