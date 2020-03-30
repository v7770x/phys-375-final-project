from constants import *
# from copy import deepcopy
from ms_functions import *
import numpy as np


'''
Add a certain amount to key values in the ODEs, useful for calculating the next y value
'''


def apply_funcs(funcs, r, y_arr):
    # rho = y_arr[0]
    # T = y_arr[1]
    # M = y_arr[2]
    # L = y_arr[3]
    # # tau = y_arr[PARAM_INDS["tau"]]

    # drho = calc_drho_dr(M,rho,r,T,L)
    # dT = calc_dT_dr(M, rho, r, T, L)
    # dM = calc_dM_dr(rho, r)
    # dL = calc_dL_dr(rho, r, T)
    # dtau = calc_dtau_dr(rho,T)
    
    
    # return np.array([drho, dT, dM, dL, dtau], float)

    return np.array([f(r, y_arr) for f in funcs], float)

'''
perform one step of rk45 integration and return the next step size sh
inputs: h = step_size, funcs = array of derivative wrapper functions, 
        y_arr  = np array of parameters w/ "rho", "T", "M", "L", "tau",
        t = indep variable"
        y = target param, tol_error = tolerance in error between 4th and 5th order approximations
outputs: h_next = next possible step size, next_y_arr = np array of updated dependent variable values
'''
def rk45_step(h, funcs, t, y_arr, tol_error):
# def rk45_step(h, funcs, params_dict, y_arr, t, tol_error):

    #####testing w/ rk4 integration########
    # k1 = h*apply_funcs(funcs, t, y_arr)
    # k2 = h*apply_funcs(funcs, t + 0.5*h, y_arr + 0.5*k1)
    # k3 = h*apply_funcs(funcs, t + 0.5*h, y_arr + 0.5*k2)
    # k4 = h*apply_funcs(funcs, t + h, y_arr + k3)
    # next_y_arr = y_arr + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0

    # return (h,next_y_arr)

    # calculate next k values
    k1 = h*apply_funcs(funcs, t,            y_arr)
    k2 = h*apply_funcs(funcs, t + 1/4*h,    y_arr + 1/4*k1)
    k3 = h*apply_funcs(funcs, t + 3/8*h,    y_arr + 3/32*k1 + 9/32*k2) 
    k4 = h*apply_funcs(funcs, t + 12/13*h,  y_arr + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3)
    k5 = h*apply_funcs(funcs, t + h,        y_arr + 439/216*k1 - 8*k2 + 3680/513*k3 -845/4104*k4)
    k6 = h*apply_funcs(funcs, t + 1/2*h,    y_arr - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40 * k5)

    # find next 4th and 5th order approximations
    next_y_arr = y_arr + 25/216*k1 + 1408/2565*k3 + 2197/4101*k4 - 1/5*k5
    next_z_arr = y_arr + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6

    # calculate minimum adaptive scaling factor
    # err = np.sqrt(np.sum(np.power(next_y_arr - next_z_arr,2)))
    # print(err)
    # h_next = h
    # if err >1:
    #     h_next = 2.0*h
    # else:
    #     h_next = h/2.0

    #find difference b/w 4th and 5th order rk approximations
    diff = np.fabs(next_z_arr - next_y_arr)

    #avoid division by 0
    is_zero = diff==0
    s = (np.fabs(next_y_arr)* tol_error/(2*(diff + is_zero )) )**(1/4) + 4*is_zero

    #make sure scaling factor does not grow or fall at an unbounded rate
    h_next = h* 0.9*max(min(np.min(s),2), 0.5)

    #make sure the step size is not too small and not too large:
    if(h_next > MAX_STEP_SIZE):
        h_next = MAX_STEP_SIZE
    if h_next < MIN_STEP_SIZE:
        h_next = MIN_STEP_SIZE
    # print("diff", diff, "s:", s, "h_next:", h_next)

    return (h_next, next_y_arr)
