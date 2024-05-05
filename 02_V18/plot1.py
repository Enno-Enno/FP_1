import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

x = np.genfromtxt("data1.txt", unpack=True)


def WQ_compton(theta):
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
