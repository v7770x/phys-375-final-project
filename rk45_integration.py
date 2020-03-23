from constants import *
import math
from copy import deepcopy


'''
Add a certain amount to key values in the ODEs, useful for calculating the next y value
'''


def add_to_vals(params, t_val, y_val):
    param_dict = deepcopy(params)
    param_dict["rho"] += y_val
    param_dict["tau"] += y_val
    param_dict["M"] += y_val
    param_dict["L"] += y_val
    param_dict["T"] += y_val
    param_dict["r"] += t_val
    return param_dict


'''
perform one step of rk45 integration and return the next step size sh
inputs: h = step_size, f = derivative function, 
        params = dictionary of parameters w/ tau, M, L, T, r, and rho,
        y = target param, tol_error = tolerance in error between 4th and 5th order approximations
outputs: h_next = next possible step size, y_next = next updated parameter value
'''
def rk45_step(h, f, params, y, tol_error):

    #####testing w/ rk4 integration########
    # k1 = h*f(params)
    # k2 = h*f(add_to_vals(params, 0.5*h, 0.5*k1))
    # k3 = h*f(add_to_vals(params, 0.5*h, 0.5*k2))
    # k4 = h*f(add_to_vals(params, h, k3))
    # y_next =  params[y] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    # return (h, y_next)

    # calculate next k values
    k1 = h*f(params)
    k2 = h*f(add_to_vals(params, 1/4*h,     1/4*k1))
    k3 = h*f(add_to_vals(params, 3/8*h,     3/32*k1 + 9/32*k2))
    k4 = h*f(add_to_vals(params, 12/13*h,   1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3))
    k5 = h*f(add_to_vals(params, h,         439/216*k1 - 8*k2 + 3680/513*k3 -845/4104*k4))
    k6 = h*f(add_to_vals(params, 1/2*h,     -8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40 * k5))

    # find next 4th and 5th order approximations
    y_next = params[y] + 25/216*k1 + 1408/2565*k3 + 2197/4101*k4 - 1/5*k5
    z_next = params[y] + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6

    # calculate adaptive scaling factor
    s = 1
    if math.fabs(z_next - y_next) > 0:
        s = (tol_error/(2*math.fabs(z_next - y_next)))**(1/4)

    #make sure scaling factor does not grow or fall at an unbounded rate
    h_next = 0.9*min(max(s, 0.3), 2)

    #make sure the step size is not too small and not too large:
    if(h_next < MAX_STEP_SIZE):
        h_next = MAX_STEP_SIZE
    if h_next < MIN_STEP_SIZE:
        h_next = MIN_STEP_SIZE

    return (h_next, y_next)

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

