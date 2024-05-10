import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import physical_constants as const
import scipy.constants as con 
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
import sympy as s
#s.init_printing()

m0csquared, unit, delta_mc  = const["electron mass energy equivalent"]
# print(m_0 *c**2, m0csquared, unit)
m0csquared_keV = m0csquared / con.kilo / con.eV

## Wichtig und richtig!

def WQ_Energie( E_gemessen,E_gamma,k, untergrund):
    # k = const * pi r_e^2 / (m_e c² epsilon²)* delta_E
    epsilon = E_gamma/m0csquared_keV
    T = E_gemessen

    t = T/E_gamma
    WQ = k * (2+ t**2/(epsilon**2 *(1-t)**2) + t /(1-t)* (t- 2/epsilon)) + untergrund
    return WQ


### ok
def E_e(E_gamma, theta):
    m_0 , unit, delta_m= const["electron mass"]
    c, unit, delta_c = const["speed of light in vacuum"]

    E_e = E_gamma *(1- 1/ (1+  E_gamma *((1-np.cos(theta))/(m0csquared_keV))) )
    return E_e



### Ignorieren
def WQ_compton(E_gamma,theta):
    ## Klein Nishima Formel 
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






def expected_N_compton(E):
    # s.diff(E_e,theta) Sympy ableitung
    diffETheta= E_gamma**2*sin(theta)/(mc*(E_gamma*(1 - cos(theta)/mc) + 1)**2)
    E_e()


def theta_inv(E_e, E_gamma):
    # m0csquared, unit, delta_mc  = const["electron mass energy equivalent"]
    return np.arccos(-(-E_gamma**2 + E_gamma*E_e + E_e*m0csquared_keV)/(E_gamma*(E_gamma - E_e))) #+ 2*np.pi 
    # np.arccos((E_gamma**2 - E_gamma*E_e - E_e*m0csquared_keV)/(E_gamma*(E_gamma - E_e)))]
    

def test_e():
    E_gamma = 1 # keV 
    print(E_e(E_gamma, 0.))
    print(E_e(E_gamma, np.pi/2))
    print(E_e(E_gamma, np.pi))
    theta = np.linspace(0,np.pi)
    plt.xlabel("theta")
    plt.ylabel("E_e in keV")
    plt.plot(theta,E_e(E_gamma,theta), label="E_e bei E_gamma = const 40 keV" )
    plt.show()
    # print(E_e(E_gamma, theta))
    # print(con.eV)
    # m0csquared, unit, delta_mc  = const["electron mass energy equivalent"]
    # m0csquared_keV = m0csquared / con.kilo/con.eV
    # print(m0csquared, m0csquared_keV)
    plt.clf()
    inverted = theta_inv(E_e(theta,E_gamma),E_gamma)
    print(np.shape(inverted),np.shape(theta))
    plt.plot(theta,inverted, label="invertierte formel")
    plt.show()

def sympy_tests():
    E_gamma, theta, m0csquared_keV,cst = s.symbols("E_gamma theta m0csquared_keV cst")
    E_e = E_gamma *(1- 1/ (1+  E_gamma *((1-s.cos(theta))/(m0csquared_keV))) )
    print(s.simplify(E_e))
    print(s.solve(E_e-cst, theta ))
    print(s.simplify(s.integrate(E_e,(theta,0,s.pi))))
    


if __name__ == "__main__":
    print("test")
    # test_e()
    # sympy_tests()



