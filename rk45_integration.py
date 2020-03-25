from constants import *
# from copy import deepcopy
import numpy as np


'''
Add a certain amount to key values in the ODEs, useful for calculating the next y value
'''


# def add_to_vals(params, t_val, y_val):
#     param_dict = deepcopy(params)
#     param_dict["rho"] += y_val
#     param_dict["tau"] += y_val
#     param_dict["M"] += y_val
#     param_dict["L"] += y_val
#     param_dict["T"] += y_val
#     param_dict["r"] += t_val
#     return param_dict


# def add(params, next_r_val, y_vals_arr):
#     param_dict = deepcopy(params)
#     for i, param in enumerate(param_dict):
#         param_dict[param] += y_vals_arr


def apply_funcs(funcs, t, y_arr):
    return np.array([f(t, y_arr) for f in funcs])

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
    # k1 = h*f(params)
    # k2 = h*f(add_to_vals(params, 0.5*h, 0.5*k1))
    # k3 = h*f(add_to_vals(params, 0.5*h, 0.5*k2))
    # k4 = h*f(add_to_vals(params, h, k3))
    # y_next =  params[y] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    # return (h, y_next)

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

# print(a)
# def test_deriv(vals):
#     return vals["M"]**2 + 1


# params = {"M": 0, "r": 0}
# step_size = 0.2
# while params["r"] < 1.41:
#     print(params, step_size)
#     next_r = params["r"]+ step_size
#     (step_size, y_next) = rk45_step(step_size, test_deriv, params, "M", 2e-5)
#     params = {"M": y_next, "r": next_r}

