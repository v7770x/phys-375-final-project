from constants import *

# this file contains the functions for the standard main sequence (the derivatives, energy generation, opacities, etc.)

###helper functions:

# P funcs
def calc_P(density,temp):
    return ((3.0*(const.pi)**2.0)**(2.0/3.0)*(const.hbar**2.0)*(density/const.m_p)**(5.0/3.0))/(5.0*const.m_e) + (density*const.k*temp)/(mu*const.m_p) + (a*(temp)**4.0)/(3.0) 

def calc_dP_drho(density,temp):
    return ((3.0*(const.pi)**2.0)**(2.0/3.0)*(const.hbar**2.0)*(density/const.m_p)**(2.0/3.0))/(3.0*const.m_e*const.m_p) + (const.k*temp)/(mu*const.m_p)

def calc_dP_dT(density,temp):
    return (density*const.k)/(mu*const.m_p) + (4.0*a*temp**3.0)/(3.0)

#kappa funcs
def kappa_ff(density,temp):
    return (1.0e24) * (1.0+X) * (Z + 0.0001)* ((density/1.0e3)**0.7) * (temp**(-3.5))

def kappa_H(density,temp):
    return (2.5e-32) * (Z / 0.02) * ((density/1.0e3)**0.5) * ((temp)**9.0)

def calc_kappa(density,temp):
    kappa = ((1.0/kappa_H(density,temp)) + (1.0/min(kappa_ff(density,temp), kappa_es)))**-1.0
    if temp > 1.0e4:
        kappa = ((1.0 / kappa_H(density, temp)) + (1.0 / max(kappa_es, kappa_ff(density,temp))))**-1.0
    return kappa

#epsilon funcs
def calc_epsilon(density,temp):
    return epsilon_PP(density,temp) + epsilon_CNO(density,temp)

def epsilon_PP(density,temp):
    return 1.07e-7*(density/1.0e5) * (X**2.0) * ((temp/1.0e6)**4.0)

def epsilon_CNO(density,temp):
    return 8.24e-26 * (density/1.0e5) * 0.03 * (X**2.0) * ((temp / (1.0e6))**19.9)
    
#deriv temp funcs
def dT_dr_radiative( mass, density, radius, temp, luminosity):
    return (3.0*calc_kappa(density,temp)*density*luminosity)/(16.0*const.pi*4.0*sigma*temp**3.0*radius**2.0)
    
def dT_dr_convective( mass, density, radius, temp, luminosity):
    pressure = calc_P(density, temp)
    return (1.0 - (1.0 / gamma)) * (temp / pressure) * ((const.G * mass * density) / (radius**2.0))


####main derivative functions:
def calc_drho_dr( mass, density, radius, temp, luminosity):
    return -((const.G*mass*density)/(radius**2.0)+calc_dP_dT(density,temp)*calc_dT_dr(mass,density,radius,temp,luminosity))/calc_dP_drho(density,temp)

def calc_dT_dr( mass, density, radius, temp, luminosity):
    return -min(abs(dT_dr_radiative(mass, density, radius, temp, luminosity)),abs(dT_dr_convective(mass,density,radius,temp,luminosity)))

def calc_dM_dr( density, radius):
    return 4.0 * const.pi * (radius**2.0) * density

def calc_dL_dr(density,radius,temp):
    return 4.0 * const.pi * (radius**2.0) * density * calc_epsilon(density,temp)

def calc_dtau_dr(density,temp):
    return calc_kappa(density,temp) * density


'''
    Wrapper functions used for generic integration function
    Pass in a dictionary of all parameters in each function, and it 
    will call the appropriate function with the appropriate parameters
 '''


def calc_dtau_dr_wrapper(r_val, param_vals):
    return calc_dtau_dr(param_vals[PARAM_INDS["rho"]], param_vals[PARAM_INDS["T"]])

def calc_dT_dr_wrapper(r_val, param_vals):
    return calc_dT_dr(param_vals[PARAM_INDS["M"]], param_vals[PARAM_INDS["rho"]], r_val, param_vals[PARAM_INDS["T"]], param_vals[PARAM_INDS["L"]])

def calc_dL_dr_wrapper(r_val, param_vals):
    return calc_dL_dr(param_vals[PARAM_INDS["rho"]], r_val, param_vals[PARAM_INDS["T"]])

def calc_dM_dr_wrapper(r_val, param_vals):
    return calc_dM_dr(param_vals[PARAM_INDS["rho"]], r_val)

def calc_drho_dr_wrapper(r_val, param_vals):
    return calc_drho_dr(param_vals[PARAM_INDS["M"]], param_vals[PARAM_INDS["rho"]], r_val, param_vals[PARAM_INDS["T"]], param_vals[PARAM_INDS["L"]])

