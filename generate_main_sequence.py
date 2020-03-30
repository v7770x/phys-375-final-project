from ms_functions import *
from constants import *
from rk45_integration import *
import matplotlib.pyplot as plt
import numpy as np

# for latex labeling
from matplotlib import rc
rc('text', usetex=True)

d_dr_functions_dict = {"rho": calc_drho_dr_wrapper, "T": calc_dT_dr_wrapper, "M": calc_dM_dr_wrapper,
                    "L": calc_dL_dr_wrapper, "tau": calc_dtau_dr_wrapper}

'''helper functions'''


# checks if tau infinity has been reached
def reached_tau_infinity(params):
    drho_dr = calc_drho_dr(params["rho"], params["T"],
                    params["r"], params["M"],params["L"])
    kappa = calc_kappa(params["rho"], params["T"])
    # if drho_dr != 0:
    #     print("delta_tau: ", kappa * params["rho"]**2/np.fabs(drho_dr), "  drho_dr: ", drho_dr, "  mult_top: ",  kappa * params["rho"]**2)
    if drho_dr != 0 and  kappa * params["rho"]**2/np.fabs(drho_dr) < tau_infinity_margin:
        return True
    return False


# function to extract the necessary parameters for the current index of integration from the dictionary of parameters
def generate_curr_dict(params, curr_index):
    return {param: params[param][curr_index] for param in params}

def generate_integrated_params_arr(params, curr_index):
    # return np.array([params["rho"][curr_index], params["T"][curr_index], 
    #     params["M"][curr_index], params["L"][curr_index], params["tau"][curr_index]],float)
    curr_arr = [0.0]*len(PARAM_INDS)
    for param in PARAM_INDS:
        curr_arr[PARAM_INDS[param]] = params[param][curr_index]
    return np.array(curr_arr, float)


# function to find the index at which tau_infinity - tau = 2/3, used to define R*
def find_r_star_index(tau_vals):
    #find tau infinity
    tau_infinity = tau_vals[-1]
    print("tau_infinity", tau_infinity)

    #find r_star index
    r_star_index = np.argmin(np.abs(tau_infinity-np.array(tau_vals) - (2.0/3.0)))
    print(r_star_index)

    if r_star_index == 0:
        return len(tau_vals) - 1
    
    return r_star_index

'''main integration function given T_c and rho_c

    returns: MS_params = dictionary of main sequence parameters solved for over range, 
            num_vals = number of values solved for'''
def solve_eqns(T_c, rho_c):
    # declare arrays/dictionaries and other necessary variables
    MS_params = {"r": [R_0], "rho": [rho_c], "T": [T_c], "M": [4*const.pi/3*R_0**3*rho_c],
         "L": [4*const.pi/3*R_0**3*rho_c*calc_epsilon(rho_c, T_c)], "tau": [calc_kappa(rho_c, T_c)*rho_c]}
    
    #create list of derivative functions
    d_dr_functions_arr = [calc_drho_dr]*len(PARAM_INDS)
    for param in PARAM_INDS:
        d_dr_functions_arr[PARAM_INDS[param]] = d_dr_functions_dict[param]

    # changing step size:
    step_size = 10000

    # keep track of number of values
    num_vals = 1 

    # variable to keep track of the current parameters
    curr_params_dict = generate_curr_dict(MS_params, num_vals - 1)

    # main integration loop, stop if mass too great or reached tau_infinity
    while curr_params_dict["M"] < 10**3 * M_sun and curr_params_dict["r"]<1.0e10 and not reached_tau_infinity(curr_params_dict):
        print(curr_params_dict, "step_size: ", step_size, "   num_vals: ", num_vals)
        # increment r
        MS_params["r"].append(curr_params_dict["r"] + step_size)

        #array of y values
        curr_params_arr = generate_integrated_params_arr(MS_params, num_vals - 1)

        #go through 1 step of rk 45 integration
        (step_size, next_params_arr) = rk45_step(step_size, d_dr_functions_arr, curr_params_dict["r"], curr_params_arr, TOL_RK_ERROR)

        #update values in dictionary
        for param in PARAM_INDS:
            MS_params[param].append(next_params_arr[PARAM_INDS[param]])

        # increment number of values, update parameter dictionary for next step of integration
        num_vals += 1
        curr_params_dict = generate_curr_dict(MS_params, num_vals - 1)
    return MS_params


# declare ICs
rho_c = 58560.0
T_c = 8.23e6

# solve the equations with the ICs
MS_params = solve_eqns(T_c, rho_c)

# find R*
r_star_index = find_r_star_index(MS_params["tau"])
r_star_params = generate_curr_dict(MS_params,  r_star_index)


#plot
for param in MS_params:
    MS_params[param] = np.array(MS_params[param][0:r_star_index])

plt.plot(MS_params["r"], MS_params["M"]/r_star_params["M"], label ="M")
plt.plot(MS_params["r"], MS_params["L"]/r_star_params["L"], label ="L")
plt.plot(MS_params["r"], MS_params["T"]/T_c, label ="T")
# plt.plot(MS_params["r"], MS_params["L"]/r_star_params["L"], label ="L")
plt.legend()
plt.show()


r_star_params["M"]/=M_sun
r_star_params["L"]/=L_sun
r_star_params["r"]/=R_sun
print("r_star_params:", r_star_params)

