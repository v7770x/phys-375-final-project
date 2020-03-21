from main_sequence_functions import *
from constants import *
from rk45_integration import *
import matplotlib.pyplot as plt
import numpy as np

# for latex labeling
from matplotlib import rc
rc('text', usetex=True)

'''helper functions'''
#checks if tau infinity has been reached
def reached_tau_infinity(params):
    drho_dr = calc_drho_dr(params["rho"], params["T"],
                    params["r"], params["M"],params["L"])
    kappa = calc_kappa(params["rho"], params["T"])
    # if drho_dr != 0:
    #     print("delta_tau: ", kappa * params["rho"]**2/math.fabs(drho_dr), "  drho_dr: ", drho_dr, "  mult_top: ",  kappa * params["rho"]**2)
    if drho_dr != 0 and  kappa * params["rho"]**2/math.fabs(drho_dr) < tau_infinity_margin:
        return True
    return False

#function to extract the necessary parameters for the current index of integration from the dictionary of parameters
def generate_curr_dict(params, wanted_params, curr_index):
    curr_dict = {}
    for param in wanted_params:
        curr_dict[param] = params[param][curr_index]
    return curr_dict

#function to find the index at which tau_infinity - tau = 2/3, used to define R*
def find_r_star_index(tau_vals):
    n_vals = len(tau_vals)
    r_star_index = n_vals - 1
    tau_infinity = tau_vals[r_star_index]

    #keep checking previous vals of tau until tau_inf - tau approx = 2/3
    while r_star_index>=0 and math.fabs(tau_infinity - tau_vals[r_star_index] - 2/3)> tau_inf_minus_tau_margin:
        r_star_index -= 1
    return r_star_index

'''main integration function given T_c and rho_c
    returns: MS_params = dictionary of main sequence parameters solved for over range, 
            num_vals = number of values solved for'''
def solve_eqns(T_c, rho_c):
    #declare arrays/dictionaries and other necessary variables
    MS_params = {"rho": [rho_c], "T": [T_c], "r": [0.0], "M": [0.0], "L":[0.0], "tau": [0.0]}
    integrated_params = ["rho", "T", "M", "L", "tau"]
    d_dr_functions = {"rho": calc_drho_dr_wrapper, "T": calc_dT_dr_wrapper, "M": calc_dM_dr_wrapper,
                    "L":calc_dL_dr_wrapper, "tau": calc_dtau_dr_wrapper}
    dict_params = ["r", "rho", "T", "M", "L", "tau"]

    #changing step size:
    step_size = 1000

    #keep track of number of values
    num_vals = 1 

    #variable to keep track of the current parameters
    curr_params = generate_curr_dict(MS_params, dict_params, num_vals - 1)
    

    #main integration loop, stop if mass too great or reched tau_infinity
    while curr_params["M"] < 10**3 * M_sun and not reached_tau_infinity(curr_params):

        # print(curr_params)
        #increment r
        MS_params["r"].append(MS_params["r"][num_vals-1] + step_size)
        
        #use minimum step size from RK 45 for next step size
        next_step_size = 1000

        #perform 1 step of integration for each parameter
        for param in integrated_params:
            # print(curr_params)
            #put through rk45
            (next_step_size_candidate, next_param_val ) = rk45_step(step_size, d_dr_functions[param], curr_params, param, tol_rk_error)
            #check for step size
            if next_step_size_candidate <next_step_size:
                next_step_size = next_step_size_candidate
            #append next calculated value
            MS_params[param].append(next_param_val)
        #increment number of values, update_step_size update parameter dictionary for next step of integration
        num_vals+=1
        step_size = next_step_size
        curr_params = generate_curr_dict(MS_params, dict_params, num_vals - 1)
    return MS_params


#declare ICs
rho_c = 58560.0
T_c = 8.23e6

#solve the equations with the ICs
MS_params = solve_eqns(T_c, rho_c)

#find R*
r_star_index = find_r_star_index(MS_params["tau"])
r_star_params  = generate_curr_dict(MS_params, ["r", "rho", "T", "M", "L", "tau"], r_star_index)
r_star_params["M"]/=M_sun
r_star_params["L"]/=L_sun
r_star_params["r"]/=R_sun
print("r_star_params:", r_star_params)

