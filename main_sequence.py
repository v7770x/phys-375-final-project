# this file consttains the calculations for the standard main sequence.

from constants import *
import numpy as np
# Scipy constants are imported in constants.py; call them with const.xxx

"""
Equation definitions:
"""


# pressure - P(rho, T)
def P(rho, T):
    return ((3*np.pi**2)**(2/3)/5)*(const.hbar**2/const.m_e)*(rho/const.m_p)**(5/3) + rho*(const.k*T/mu/const.m_p) + (1/3)*a*T**4


# dP/drho
def dP_drho(rho, T):
    return ((3*np.pi**2)**(2/3)/3)*(const.hbar**2/const.m_e/const.m_p)*(rho/const.m_p)**(2/3) + const.k*T/mu/const.m_p


# dP/dT
def dP_dT(rho, T):
    return (rho*const.k/mu/const.m_p)+(4/3)*a*T**3


# PP-chain energy
def e_pp(rho, T):
    return 1.07e-7*(rho/1e5)*X**2*(T/1e6)**4


# CNO cycle energy
def e_cno(rho, T):
    return 8.24e-26*(rho/1e5)*X*Xcno*(T/1e6)**19.9


# total solar energy generation
def epsilon(rho, T):
    return e_pp(rho, T) + e_cno(rho, T)


# ff opacity (k_ff)
def kappa_ff(rho, T):
    return (1e24)*(Z + 0.0001)*(rho/1e3)**0.7*T**-3.5


# H- opacity (k_H)
def kappa_H(rho, T):
    return (2.5e-32)*(Z/0.02)*(rho/1e3)**0.5*T**9


# overall opacity
def kappa(rho, T):
    return ((1/kappa_H(rho, T))+(1/max(kappa_es, kappa_ff(rho, T))))**-1


# r derivative of Tau
def tau_r(rho, Kap):
    return Kap*rho


# r derivative of Luminousity
def L_r(rho, T, r, e):
    return 4*np.pi*r**2*rho*e


# r derivative of mass
def M_r(rho, r):
    return 4*np.pi*r**2*rho


# r derivative of rho
def rho_r(rho, T, r, M, dP_drho, dP_dT, dT_dr):
    return -(const.G*M*rho/r**2 + dP_dT*dT_dr)/dP_drho


# r derivative of T
def dT_dr(rho, T, r, M, L, P, Kap):
    return - min(3*Kap*rho*L/16/np.pi/a/const.c/T**3/r**2, (1 - 1/gamma)*(T/P)*const.G*M*rho/r**2)

