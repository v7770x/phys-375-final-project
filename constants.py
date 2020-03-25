# all of the constants to be used in the project
# currently using scipy.constants for most, but our own constants (X etc.) should be defined here

import scipy.constants as const
import enum

M_sun = 1.989e30
R_sun = 696340000
L_sun = 3.828e26
sigma = 5.67e-8
a = 4*sigma/const.c
X = 0.734  # the X, Y, Z constants were taken from the OPAL database, and are placeholders until we can find better.
Y = 0.25
Z = 0.016
gamma = 5/3
X_CNO = 0.03*X
mu = (2*X + 0.75*Y + 0.5*Z)**-1
kappa_es = 0.02*(1 + X)

#error tolerance between 4th and 5th order rk
TOL_RK_ERROR = 1e-6

#other error tolerances
tau_infinity_margin = 1e-1
tau_inf_minus_tau_margin = 1e-3

#step size bounds:
MIN_STEP_SIZE = 1000 #m
MAX_STEP_SIZE = 1e6 #m

#keep track of parameter indicies
PARAM_INDS = {"rho": 0, "T": 1, "M": 2, "L": 3, "tau": 4}
"""
Scipy constants: https://docs.scipy.org/doc/scipy/reference/constants.html

use those if you can, eg. const.pi for pi, or const.hbar for hbar.
"""


