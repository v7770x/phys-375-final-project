# this file contains the functions for the standard main sequence (the derivatives, energy generation, opacities, etc.)

from constants import *
# Scipy constants are imported in constants.py; call them with const.xxx

"""
Equation definitions:
"""

# pressure - P(rho, T)
def calc_P(rho, T):
    return ((3.0*(const.pi)**2.0)**(2.0/3.0)/5.0)*(const.hbar**2.0/const.m_e)*(rho/const.m_p)**(5.0/3.0) + rho*(const.k*T/(mu*const.m_p)) + (1.0/3.0)*a*(T**4.0)


# dP/drho
def calc_dP_drho(rho, T):
    return ((3.0*(const.pi)**2.0)**(2.0/3.0)/3.0)*(const.hbar**2/(const.m_e/const.m_p))*(rho/const.m_p)**(2/3) + const.k*T/(mu*const.m_p)


# dP/dT
def calc_dP_dT(rho, T):
    return (rho*const.k)/(mu*const.m_p) + (4.0/3.0)*a*(T**3.0)


# PP-chain energy
def calc_epsilon_pp(rho, T):
    return 1.07e-7*(rho/1e5)*(X**2.0)*((T/1.0e6)**4.0)


# CNO cycle energy
def calc_epsilon_cno(rho, T):
    return 8.24e-26*(rho/1e5)*X*X_CNO*(T/1e6)**19.9


# total solar energy generation
def calc_epsilon(rho, T):
    return calc_epsilon_pp(rho, T) + calc_epsilon_cno(rho, T)


# ff opacity (k_ff)
def calc_kappa_ff(rho, T):
    return (1e24)*(Z + 0.0001)* ((rho/1e3)**0.7) * (T**(-3.5))


# H- opacity (k_H)
def calc_kappa_H(rho, T):
    return (2.5e-32)*(Z/0.02)* ((rho/1.0e3)**0.5)* (T**9.0)


# overall opacity
def calc_kappa(rho, T):
    if T > 1.0e4:
        return ((1.0/calc_kappa_H(rho, T))+(1.0/max(kappa_es, calc_kappa_ff(rho, T))))**-1.0
    else:
        return ((1.0/calc_kappa_H(rho, T))+(1.0/min(kappa_es, calc_kappa_ff(rho, T))))**-1.0


# r derivative of Tau
def calc_dtau_dr(rho, T):
    return calc_kappa(rho, T)*rho


# r derivative of T
def calc_dT_dr(rho, T, r, M, L):
    if r == 0:
        return 0.0
    dT_dr_convective = (1.0 - (1.0 / gamma)) * (T / calc_P(rho, T)) * ((const.G * M * rho) / (r**2.0))
    dT_dr_radiative = (3.0*calc_kappa(rho, T)*rho*L)/(16.0*const.pi*4.0*sigma*T**3.0*r**2.0)
    return -min(abs(dT_dr_convective), abs(dT_dr_radiative))


# r derivative of Luminousity
def calc_dL_dr(rho, T, r):
    return 4.0*const.pi*(r**2)*rho*(calc_epsilon(rho, T))


# r derivative of mass
def calc_dM_dr(rho, r):
    return 4.0*const.pi*(r**2.0)*rho


# r derivative of rho
def calc_drho_dr(rho, T, r, M, L):
    if r==0: 
        return 0.0
    return -( (const.G*M*rho)/(r**2.0) + calc_dP_dT(rho, T)*calc_dT_dr(rho, T, r, M, L))/calc_dP_drho(rho,T)


'''
    Wrapper functions used for generic integration function
    Pass in a dictionary of all parameters in each function, and it 
    will call the appropriate function with the appropriate parameters
 '''


def calc_dtau_dr_wrapper(r_val, param_vals):
    return calc_dtau_dr(param_vals[PARAM_INDS["rho"]], param_vals[PARAM_INDS["T"]])

def calc_dT_dr_wrapper(r_val, param_vals):
    return calc_dT_dr(param_vals[PARAM_INDS["rho"]], param_vals[PARAM_INDS["T"]], r_val, param_vals[PARAM_INDS["M"]], param_vals[PARAM_INDS["L"]])

def calc_dL_dr_wrapper(r_val, param_vals):
    return calc_dL_dr(param_vals[PARAM_INDS["rho"]], param_vals[PARAM_INDS["T"]], r_val)

def calc_dM_dr_wrapper(r_val, param_vals):
    return calc_dM_dr(param_vals[PARAM_INDS["rho"]], r_val)

def calc_drho_dr_wrapper(r_val, param_vals):
    return calc_drho_dr(param_vals[PARAM_INDS["rho"]], param_vals[PARAM_INDS["T"]], r_val, param_vals[PARAM_INDS["M"]], param_vals[PARAM_INDS["L"]])