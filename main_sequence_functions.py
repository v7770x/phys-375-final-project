# this file contains the functions for the standard main sequence (the derivatives, energy generation, opacities, etc.)

from constants import *
# Scipy constants are imported in constants.py; call them with const.xxx

"""
Equation definitions:
"""

# pressure - P(rho, T)
def calc_P(rho, T):
    return ((3*const.pi**2)**(2/3)/5)*(const.hbar**2/const.m_e)*(rho/const.m_p)**(5/3) + rho*(const.k*T/mu/const.m_p) + (1/3)*a*T**4


# dP/drho
def calc_dP_drho(rho, T):
    return ((3*const.pi**2)**(2/3)/3)*(const.hbar**2/const.m_e/const.m_p)*(rho/const.m_p)**(2/3) + const.k*T/mu/const.m_p


# dP/dT
def calc_dP_dT(rho, T):
    return (rho*const.k/mu/const.m_p)+(4/3)*a*T**3


# PP-chain energy
def calc_epsilon_pp(rho, T):
    return 1.07e-7*(rho/1e5)*X**2*(T/1e6)**4


# CNO cycle energy
def calc_epsilon_cno(rho, T):
    return 8.24e-26*(rho/1e5)*X*X_CNO*(T/1e6)**19.9


# total solar energy generation
def calc_epsilon(rho, T):
    return calc_epsilon_pp(rho, T) + calc_epsilon_cno(rho, T)


# ff opacity (k_ff)
def calc_kappa_ff(rho, T):
    return (1e24)*(Z + 0.0001)*(rho/1e3)**0.7*T**-3.5


# H- opacity (k_H)
def calc_kappa_H(rho, T):
    return (2.5e-32)*(Z/0.02)*(rho/1e3)**0.5*T**9


# overall opacity
def calc_kappa(rho, T):
    return ((1/calc_kappa_H(rho, T))+(1/max(kappa_es, calc_kappa_ff(rho, T))))**-1


# r derivative of Tau
def calc_dtau_dr(rho, T):
    return calc_kappa(rho,T)*rho


# r derivative of T
def calc_dT_dr(rho, T, r, M, L):
    if r == 0:
        return 0.0
    return - min(3*calc_kappa(rho,T)*rho*L/16/const.pi/a/const.c/T**3/(r**2), (1 - 1/gamma)*(T/calc_P(rho,T)*const.G*M*rho/(r**2)))


# r derivative of Luminousity
def calc_dL_dr(rho, T, r):
    return 4*const.pi*(r**2)*rho*(calc_epsilon(rho,T))


# r derivative of mass
def calc_dM_dr(rho, r):
    return 4*const.pi*(r**2)*rho


# r derivative of rho
def calc_drho_dr(rho, T, r, M, L):
    if r==0: 
        return 0.0
    return -(const.G*M*rho/r**2 + calc_dP_dT(rho, T)*calc_dT_dr(rho, T, r, M, L))/calc_dP_drho(rho,T)


'''Wrapper functions used for generic integration function
    Pass in a dictionary of all parameters in each function, and it 
    will call the appropriate function with the appropriate parameters
 '''
def calc_dtau_dr_wrapper(params):
    return calc_dtau_dr(params["rho"], params["T"])

def calc_dT_dr_wrapper(params):
    return calc_dT_dr(params["rho"], params["T"], params["r"], params["M"], params["L"])

def calc_dL_dr_wrapper(params):
    return calc_dL_dr(params["rho"], params["T"], params["r"])

def calc_dM_dr_wrapper(params):
    return calc_dM_dr(params["rho"], params["r"])

def calc_drho_dr_wrapper(params):
    return calc_drho_dr(params["rho"], params["T"], params["r"], params["M"], params["L"])