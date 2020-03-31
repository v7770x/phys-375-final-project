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
    if np.isnan(drho_dr) or (drho_dr != 0 and  kappa * params["rho"]**2/np.fabs(drho_dr) < tau_infinity_margin):
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
    tau_inf_index = len(tau_vals)-1
    if np.isnan(tau_vals[tau_inf_index]):
        tau_inf_index = len(tau_vals) -2
    tau_infinity = tau_vals[tau_inf_index]
    # print("tau_infinity", tau_infinity)

    #find r_star index
    # print(np.min(np.abs(tau_infinity-np.array(tau_vals) - (2.0/3.0))))
    r_star_index = np.argmin(np.abs(tau_infinity-np.array(tau_vals[0:tau_inf_index]) - (2.0/3.0)))
    # print(r_star_index)

    if r_star_index == 0:
        return tau_inf_index
    
    return r_star_index

def get_r_star_params(MS_params):
    r_star_index = find_r_star_index(MS_params["tau"])
    r_star_params = generate_curr_dict(MS_params,  r_star_index)
    return (r_star_params, r_star_index)


def update_non_integrated_params(MS_params, next_params_arr):
    T = next_params_arr[PARAM_INDS["T"]]
    M = next_params_arr[PARAM_INDS["M"]]
    rho = next_params_arr[PARAM_INDS["rho"]]
    L = next_params_arr[PARAM_INDS["L"]]
    r = MS_params["r"][-1]

    MS_params["kappa"].append(calc_kappa(rho,T))
    MS_params["P"].append(calc_P(rho,T))
    MS_params["dL_dr"].append(calc_dL_dr(rho,r, T))
        
def f_rho_c(MS_params):
    #find surface params
    (r_star_params, r_star_index) = get_r_star_params(MS_params)

    #find f
    numerator = (r_star_params["L"] - 4.0*const.pi*sigma*(r_star_params["r"]**2.0)*(r_star_params["T"]**4.0))
    denominator = np.sqrt(4.0*const.pi*sigma*(r_star_params["r"]**2.0)*(r_star_params["T"]**4.0)*r_star_params["L"])
    return numerator/denominator

'''main integration function given T_c and rho_c

    returns: MS_params = dictionary of main sequence parameters solved for over range, 
            num_vals = number of values solved for'''
def solve_eqns(T_c, rho_c):
    # declare arrays/dictionaries and other necessary variables
    MS_params = {"r": [R_0], "rho": [rho_c], "T": [T_c], "M": [4*const.pi/3*R_0**3*rho_c],
         "L": [4*const.pi/3*R_0**3*rho_c*calc_epsilon(rho_c, T_c)], "tau": [calc_kappa(rho_c, T_c)*rho_c]
         , "kappa": [calc_kappa(rho_c, T_c)], "P": [calc_P(rho_c,T_c)], "dL_dr": [calc_dL_dr(rho_c, R_0, T_c)]}
    
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
        # print(curr_params_dict, "step_size: ", step_size, "   num_vals: ", num_vals)
        # increment r
        MS_params["r"].append(curr_params_dict["r"] + step_size)

        #array of y values
        curr_params_arr = generate_integrated_params_arr(MS_params, num_vals - 1)

        #go through 1 step of rk 45 integration
        (step_size, next_params_arr) = rk45_step(step_size, d_dr_functions_arr, curr_params_dict["r"], curr_params_arr, TOL_RK_ERROR, T_c)

        #update integrated values in dictionary
        for param in PARAM_INDS:
            MS_params[param].append(next_params_arr[PARAM_INDS[param]])
        
        #update all values in dictionary
        update_non_integrated_params(MS_params, next_params_arr)
        
        # increment number of values, update parameter dictionary for next step of integration
        num_vals += 1
        curr_params_dict = generate_curr_dict(MS_params, num_vals - 1)

    #find the surface parameters, and clip the MS_param vals to the surface, make into np arrs
    (r_star_params, r_star_index) = get_r_star_params(MS_params)
    # print(r_star_params)
    # print(MS_params["tau"])
    for param in MS_params:
        MS_params[param] = np.array(MS_params[param][0:r_star_index])

    return MS_params

''' 
Bisection to find rho_c given T_c
'''


#testing function
def f_rho_behavior():
    T_c = 8.23e6
    rho_list = []
    f_rho_list = []
    for rho_c in range(300, 500000, 2000):
        MS_params = solve_eqns(T_c, rho_c)
        r_star_index = find_r_star_index(MS_params["tau"])
        r_star_params = generate_curr_dict(MS_params,  r_star_index)
        f_rho = f_rho_c(MS_params)
        print("rho_c:",rho_c, "F_rhoc = ", f_rho)
        if not np.isnan(f_rho):
            rho_list.append(rho_c)
            f_rho_list.append(f_rho)
    plt.plot(np.array(rho_list), np.array(f_rho_list), "-b")
    plt.show()
    return (rho_list, f_rho_list)

# (rho_list, f_rho_list) = f_rho_behavior()

#given a T_c find rho_c and the star surface parameters as well as all the data for the star
#over the generation of the star
def find_rho_c_params(T_c):
    #find params at minimum and maximum rho_c, midpoint
    rho_c_lb_params = solve_eqns(T_c, RHO_C_MIN) 
    rho_c_ub_params = solve_eqns(T_c, RHO_C_MAX)
    rho_c_mid_params = solve_eqns(T_c, (RHO_C_MAX + RHO_C_MIN)/2)

    #rho tolerances
    RHO_C_DIFF_TOL = 1e-2
    MAX_NUM_BISECTIONS = 40
    RHO_C_F_TOL = 1e-2

    #keep track of num bisections
    num_bisections = 0

    #calculate f_rho_c at midpoint
    f_mid = f_rho_c(rho_c_mid_params)

    #loop through until rho_c not changing much, f_rho_c less than tol, or num_bisections too large
    while (abs(rho_c_lb_params["rho"][0] - rho_c_ub_params["rho"][0])>RHO_C_DIFF_TOL 
        and abs(f_mid)>RHO_C_F_TOL and num_bisections<MAX_NUM_BISECTIONS):

        print("num_bisections: ", num_bisections, "rho_c:", rho_c_mid_params["rho"][0], "f_rho_c", f_mid)

        if np.isnan(f_mid) or f_mid > 0:
            rho_c_ub_params = rho_c_mid_params
        elif f_mid<0:
            rho_c_lb_params = rho_c_mid_params
        # else:
        #     break
        rho_c_mid_params = solve_eqns(T_c, (rho_c_ub_params["rho"][0] + rho_c_lb_params["rho"][0])/2)

        #increment num bisections
        num_bisections +=1

        #update f_rho_c at midpoint
        f_mid = f_rho_c(rho_c_mid_params)

        # print("diff", abs(rho_c_lb_params["rho"][0] - rho_c_ub_params["rho"][0]), "f_mid", f_mid)

    #ensure the smallest f_mid value is chosen
    f_ub = abs(f_rho_c(rho_c_ub_params))
    f_lb = abs(f_rho_c(rho_c_lb_params))
    f_mid = abs(f_mid)
    if f_mid > RHO_C_F_TOL and f_ub < f_mid and f_ub < f_lb:
        rho_c_mid_params = rho_c_ub_params
        print("ub")
    elif f_mid > RHO_C_F_TOL and f_lb < f_mid and f_lb < f_ub:
        rho_c_mid_params = rho_c_lb_params
        print("lb")

    #find star params and print out
    f_mid_final = f_rho_c(rho_c_mid_params)
    rho_c_found = rho_c_mid_params["rho"][0]
    print("rho_c_final = ", rho_c_found, "T_c", rho_c_mid_params["T"][0], "f_rho_c", f_mid_final)
    (r_star_params, r_star_index) = get_r_star_params(rho_c_mid_params)
    print("generated star params: ", r_star_params)

    #print warning if bisection limit hit
    if num_bisections == MAX_NUM_BISECTIONS:
        print("WARNING: max bisection limit hit")

    #print warning if f_rho_mid no good:
    if abs(f_mid_final>10):
        print("WARNING: bisection function may not have found root, f_rho_c too high")

    #print warning if rho_c near limits 
    if rho_c_found > RHO_C_MAX - 1.0:
        print("WARNING: RHO_C TOO HIGH, MAY NOT BE REAL STAR")
    elif rho_c_found<RHO_C_MIN + 1.0:
        print("WARNING: RHO_C TOO LOW, MAY NOT BE REAL STAR")
    

    #return
    return rho_c_mid_params


'''Plotting functions'''
def generate_plots(MS_params):
    (r_star_params, r_star_index) = get_r_star_params(MS_params)
    T_c = MS_params["T"][0]
    rho_c = MS_params["rho"][0]

    #plot r vs M, L, rho, T
    plt.plot(MS_params["r"], MS_params["M"]/r_star_params["M"], label ="M")
    plt.plot(MS_params["r"], MS_params["L"]/r_star_params["L"], label ="L")
    plt.plot(MS_params["r"], MS_params["T"]/T_c, label ="T")
    plt.plot(MS_params["r"], MS_params["rho"]/rho_c, label ="rho")
    plt.legend()

    #plot r vs P
    plt.figure(2)
    plt.plot(MS_params["r"], MS_params["P"]/ MS_params["P"][0], label= "P")
    plt.legend()

    #plot r vs kappa
    plt.figure(3)
    plt.plot(MS_params["r"], np.log10(MS_params["kappa"]), label= "kappa")
    plt.legend()

    plt.show()

rho_c_params = find_rho_c_params(8.23e6)
generate_plots(rho_c_params)
# generate_plots(solve_eqns(8.23e6, 58000.0))



