import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

polarization_angles_deg = np.array(
    [
        0,
        10,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        70,
        80,
        90,
        100,
        110,
        120,
        130,
        140,
        150,
        160,
        170,
        180,
    ]
)
Emin = np.array(
    [
        [1.35, 1.36, 1.36],
        [0.94, 0.95, 0.96],
        [0.46, 0.45, 0.45],
        [0.29, 0.30, 0.31],
        [0.21, 0.20, 0.24],
        [0.17, 0.18, 0.18],
        [0.12, 0.12, 0.12],
        [0.10, 0.10, 0.10],
        [0.10, 0.10, 0.11],
        [0.12, 0.12, 0.12],
        [0.18, 0.18, 0.18],
        [0.42, 0.41, 0.42],
        [0.92, 0.92, 0.92],
        [1.60, 1.62, 1.62],
        [1.67, 1.70, 1.68],
        [1.06, 1.07, 1.10],
        [0.69, 0.63, 0.65],
        [0.36, 0.40, 0.39],
        [0.42, 0.44, 0.43],
        [0.59, 0.67, 0.66],
        [0.85, 0.86, 0.89],
        [1.13, 1.18, 1.18],
        [1.38, 1.39, 1.39],
    ]
)
Emax = np.array(
    [
        [1.57, 1.61, 1.62],
        [1.32, 1.33, 1.33],
        [1.44, 1.46, 1.42],
        [1.30, 1.34, 1.36],
        [1.25, 1.18, 1.18],
        [1.16, 1.16, 1.16],
        [1.12, 1.17, 1.16],
        [1.30, 1.31, 1.26],
        [1.36, 1.36, 1.37],
        [1.50, 1.51, 1.49],
        [1.61, 1.63, 1.61],
        [2.04, 2.00, 1.99],
        [2.12, 2.16, 2.18],
        [2.21, 2.20, 2.19],
        [2.64, 2.69, 2.69],
        [3.92, 4.00, 3.89],
        [4.73, 4.54, 4.77],
        [4.86, 4.89, 5.07],
        [5.08, 5.19, 5.07],
        [4.27, 4.16, 4.40],
        [3.60, 3.28, 3.29],
        [2.22, 2.42, 2.40],
        [1.67, 1.66, 1.66],
    ]
)

RefIndexGlasNumbers = np.array([35, 31, 30, 37, 38, 37, 37, 34, 33, 37, 34])
pressures_mbar = np.zeros(21)
pi = 0
for i, _ in enumerate(pressures_mbar):
    pressures_mbar[i] = pi
    pi += 50

counter_air = np.array(
    [
        [0, 2, 4, 6, 8, 10, 15, 17, 19, 21, 22, 25, 27, 29, 31, 34, 36, 38, 40, 42, 44],
        [0, 2, 4, 6, 8, 11, 13, 15, 17, 19, 21, 23, 25, 28, 30, 33, 35, 37, 39, 41, 43],
        [0, 2, 4, 6, 9, 10, 13, 15, 17, 19, 21, 23, 25, 27, 30, 32, 34, 36, 38, 40, 42],
        [0, 2, 4, 6, 8, 11, 13, 15, 17, 19, 21, 23, 25, 27, 30, 32, 34, 36, 38, 40, 42],
    ]
)

T_cel = 20.1

if __name__ == "__main__":
    # print("shape(RefIndexGlasNumbers)", np.shape(RefIndexGlasNumbers))
    # print(
        # "shape(Emax),shape(Emin),shape(polarization_angles_deg)",
        # np.shape(Emax),
        # np.shape(Emin),
        # np.shape(polarization_angles_deg),
    # )
    # print(
        # "np.shape(counter_air),np.shape(pressures_mbar)",
        # np.shape(counter_air),
        # np.shape(pressures_mbar),
    # )
    for i, deg in enumerate(polarization_angles_deg):
        print(deg," & ", Emin[i,0]," & ", Emax[i,0], r" \\")
        print(" & ", Emin[i,1]," & ", Emax[i,1], r" \\")
        print(" & ", Emin[i,2]," & ", Emax[i,2], r" \\")
