import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
import data as d
import pprint

pp = pprint.PrettyPrinter(indent=4)


# shape(RefIndexGlasNumbers) (11,)
# shape(Emax),shape(Emin),shape(polarization_angles_deg) (23, 3) (23, 3) (23,)
# np.shape(counter_air),np.shape(pressures_mbar) (4, 21) (21,)

### Interferometer optimaler Polarisationswinkel
angles = d.polarization_angles_deg

EminDelta = np.zeros(len(angles))
EminNom = np.zeros(len(angles))
EmaxDelta = np.zeros(len(angles))
EmaxNom = np.zeros(len(angles))
for i, _ in enumerate(angles):
    EminNom[i] = np.mean(d.Emin[i, :])
    EminDelta_i = np.std(d.Emin[i, :])
    if EminDelta_i < 0.005:
        EminDelta[i] = 0.005
    else:
        EminDelta[i] = EminDelta_i
    
    EmaxNom[i] = np.mean(d.Emax[i, :])
    EmaxDelta_i = np.std(d.Emax[i, :])
    if EmaxDelta_i < 0.005:
        EmaxDelta[i] = 0.005
    else:
        EmaxDelta[i] = EmaxDelta_i
Emax = unp.uarray(EmaxNom, EmaxDelta)
Emin = unp.uarray(EminNom, EminDelta)

contrast= lambda Emin,Emax: (Emax-Emin)/(Emax+Emin)

consNom = np.zeros(len(Emin))
consDel = np.zeros(len(Emin))
for i,_ in enumerate(Emin):
    con = contrast(Emin[i],Emax[i])
    consNom[i] = unc.nominal_value(con)
    consDel[i] = unc.std_dev(con)
cons = unp.uarray(consNom,consDel)


p = False
if p:
    for i,_ in enumerate(Emin):
        print(angles[i]," & ",Emax[i]," & ", Emin[i], cons[i]," & ", cons[i]," \\\\")

degReturnRad = lambda deg: 2 *np.pi/360 *deg

def expected_contrast_deg(theta):
    theta_rad = degReturnRad(theta)
    return np.abs(2* np.cos(theta_rad)* np.sin(theta_rad))

theta_x = np.linspace(0,180,180)
plt.plot(theta_x, expected_contrast_deg(theta_x), label="Expected contrast")
plt.errorbar(angles, consNom,fmt="x", yerr=consDel, label="Measured contrasts")
plt.grid()
plt.legend()
plt.xlabel("theta/°")
# plt.show()

plt.savefig("build/plot1.pdf")