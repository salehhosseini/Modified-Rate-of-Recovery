''' This code is to derive the analytical rate of recovery based on powered
simplest walking model as representing the portion fo the energy wasted
during the step-to-step transition'''

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize      import fsolve

# velocity range
v_ave    = np.linspace(0.8, 1.6, 100)

###
## initial guess
alpha_ini = 0.5 * v_ave ** 0.42

## finding alpha when v_ave = 1.5m/s
def find_index_larger_than(arr, target):
    for i in range(len(arr)):
        if arr[i] > target:
            return i
    # If no element in the array is larger than the target
    return -1

v_target = find_index_larger_than(v_ave, 1.5)
c = 1 / (alpha_ini[v_target] / 0.4)

## Adjusted Alpha
alpha = 0.5 * c * v_ave**0.42
x_ticks = [0.8, 1.0, 1.2, 1.4, 1.6]

plt.figure(1, figsize = (14, 6))
plt.subplot(1,3,1)
plt.plot(v_ave, alpha * 2, color = 'k', linewidth = 2.5)

plt.xlabel('Walking Velocity', fontsize = 12.5); plt.xticks(x_ticks, fontsize = 12.5)
plt.ylabel('m', fontsize = 15); plt.yticks(fontsize = 12.5)
plt.title('A. Step Length', loc = 'left', fontsize = 15)

##
E_total  = 0.5 * v_ave**2
PO       = 0.125 * c**2 * v_ave**2.84
Energy_Recovery = 0.5 * v_ave**2 - 0.125 * c**2 * v_ave**2.84

plt.subplot(1,3,2); 
plt.plot(v_ave, E_total, color = 'k', linewidth = 2.5, label = 'Total Energy')
plt.plot(v_ave, PO, color = 'navy', linewidth = 2.5, label = 'Push-Off Work')
plt.plot(v_ave, Energy_Recovery, 'darkred', linewidth = 2.5, label = 'Recovery')
plt.legend(frameon = False, fontsize = 10)
plt.xlabel('Walking Velocity', fontsize = 12.5); plt.xticks(x_ticks, fontsize = 12.5)
plt.ylabel('J/kg', fontsize = 15); plt.yticks(fontsize = 12.5)
plt.title('B. Mechanical Energies', loc = 'left', fontsize = 15)

##
Recovery_rate = 1 - (np.tan(0.5 * c * v_ave**0.42))**2
Recovery_rate_app = 1 - 0.25 * c**2 * v_ave**0.84
Dissipation_rate = 0.25 * c**2 * v_ave**0.84

plt.subplot(1,3,3)
plt.plot(v_ave, Recovery_rate_app * 100, color = 'grey', linewidth = 2.5, label = 'Rate of Reacovery (App.)', linestyle = '--')
plt.plot(v_ave, Recovery_rate * 100,     color = 'darkred', linewidth = 2.5, label = 'Rate of Reacovery')
plt.plot(v_ave, Dissipation_rate * 100, color = 'navy', linewidth = 2.5, label = 'Dissipation_Rate')
plt.legend(frameon = False, fontsize = 10)
plt.xlabel('Walking Velocity', fontsize = 12.5); plt.xticks(x_ticks, fontsize = 12.5)
plt.ylabel('%', fontsize = 15); plt.yticks(fontsize = 12.5)
plt.title('C. Rates', loc = 'left', fontsize = 15)

plt.tight_layout(pad = 0.5, w_pad=0.25)

#### Rate of recovery when the push-off is suboptimal
## Initial conditions
v_ave     = 1.25                    # average walking velocity
alpha_ave = 0.5 * c * v_ave**0.42   # angle between legs

# optimal push-off impulse for a give walking velocity
po_opt = v_ave * np.tan(alpha_ave)

# Range of push-off
po_rand = np.linspace(0, po_opt, 100)

# mid-transition velocity angle with v_ave
gama = np.arctan(po_rand / v_ave)

# Collision for a delayed push-off
co = v_ave * np.sin(2 * alpha_ave) - po_rand * np.cos(2 * alpha_ave)

#  Collision work: 0.5 * co**2
CO_rand = 0.5 * co**2

# Push-off work: 0.5 * po**2
PO_rand = 0.5 * po_rand**2

## Rate of return
# RoR = 1 - (np.sin(2 * alpha_ave) - (po_rand / v_ave) * np.cos(2 * alpha_ave))** 2
RoR     = 1 - (np.sin(2 * alpha_ave - gama))**2 / (np.cos(gama))**2
RoR_app = 1 - (2 * alpha_ave - gama)**2

## 
po_fract = np.linspace(0, 100, 100)

## Plotting
plt.figure(2, figsize=(14,6))
#
plt.subplot(1,4,1)
plt.plot(po_fract, alpha_ave * np.linspace(1, 1, 100), color = 'k', linestyle = '--', label = 'Nominal Angle')
plt.plot(po_fract, gama, color = 'k', linewidth = 2.5, label = 'Random Push-off Angle')
plt.legend(frameon = False)
plt.xticks(fontsize = 12.5); plt.xlabel('% PO_nominal', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('Rad.', fontsize = 15)
plt.title('A. Mid-Transition Angle', loc = 'left', fontsize = 15)

#
plt.subplot(1,4,2)
plt.plot(po_fract, po_rand, color = 'darkred', label = 'Push-off Impulse', linewidth = 2.5)
plt.plot(po_fract, co, color = 'navy', label = 'Colliion Impulse', linewidth = 2.5)
plt.legend(frameon = False)
plt.xticks(fontsize = 12.5); plt.xlabel('% PO_nominal', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('m/s', fontsize = 15)
plt.title('B. Transition Impulses', loc = 'left', fontsize = 15)

#
plt.subplot(1,4,3)
plt.plot(po_fract, PO_rand, color = 'darkred', label = 'Push-Off Work', linewidth = 2.5)
plt.plot(po_fract, CO_rand, color = 'navy', label = 'Collision Work', linewidth = 2.5)
plt.legend(frameon = False)
plt.xticks(fontsize = 12.5); plt.xlabel('% PO_nominal', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('J/kg', fontsize = 15)
plt.title('C. Transition Works', loc = 'left', fontsize = 15)

# 
plt.subplot(1,4,4)
plt.plot(po_fract, RoR_app * 100, color = 'grey', linewidth = 2.5, label = 'Rate of Reacovery (App.)', linestyle = '--')
plt.plot(po_fract, RoR * 100, color = 'k', linewidth = 2.5, label = 'Rate of Reacovery')
plt.xticks(fontsize = 12.5); plt.xlabel('% PO_nominal', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('%', fontsize = 15)
plt.title('D. Rate of Recovery', loc = 'left', fontsize = 15)
plt.legend(frameon = False)

plt.tight_layout()

####
dh_range = np.linspace(0, 0.05, 1000)
g = 9.81

# Compute v_post_opt
v_post_opt = np.sqrt(v_ave**2 + 2 * g * dh_range)

# Define the function based on the equation V_post_opt = V_ave * cos(2alpha - gama) / cos(gama)
def equation_to_solve(gama_opt, v_post_opt, alpha_ave, v_ave):
    return v_post_opt - (v_ave * np.cos(2 * alpha_ave - gama_opt) / np.cos(gama_opt))

# Initialize a list to store the solutions
gama_opt_list = []

for v_post in v_post_opt:
    # Initial guess for gama (starting point for the solver)
    gama_guess = 0.35  # Adjust as necessary based on expected range for gama

    # Use fsolve to find the root of the equation
    gama_opt_solution = fsolve(equation_to_solve, gama_guess, args=(v_post, alpha_ave, v_ave))

    # Append the solution to the list
    gama_opt_list.append(gama_opt_solution[0])

# Convert the list to a NumPy array
gama_opt_array = np.array(gama_opt_list)

## Gama_opt must be less than 2 * alpha
gama_opt_array = gama_opt_array[gama_opt_array < 2 * alpha_ave]

RoR_uneven     = 1 - (np.sin(2 * alpha_ave - gama_opt_array))**2 / (np.cos(gama_opt_array))**2
RoR_uneven_app = 1 - (2 * alpha_ave - gama_opt_array)**2

po_up = v_ave * np.tan(gama_opt_array); PO_up = 0.5 * po_up**2
co_up = v_ave * (np.sin(2 * alpha_ave - gama_opt_array)) / (np.cos(gama_opt_array)); CO_up = 0.5 * co_up**2

##
plt.figure(3, figsize=(14, 6)) 

plt.subplot(1,3,1)
plt.plot(dh_range, gama_opt_array, linewidth = 2.5, color = 'k')
plt.title('A. Mid-Transition Angle', loc = 'left', fontsize = 15)

plt.xticks(fontsize = 12.5); plt.xlabel('Step-up Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('Rad.', fontsize = 15)

plt.subplot(1,3,2)
plt.plot(dh_range, PO_up, label = 'Optimal Push-off', linewidth = 2.5, color = 'darkred')
plt.plot(dh_range, CO_up, label = 'Optimal Collision', linewidth = 2.5, color = 'navy')

plt.title('B. Step Transition Works', loc = 'left', fontsize = 15)
plt.legend(frameon = False)
plt.xticks(fontsize = 12.5); plt.xlabel('Step-up Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('J/kg', fontsize = 15)

plt.subplot(1,3,3)
plt.plot(dh_range, RoR_uneven * 100, linewidth = 2.5, color = 'k', label = 'Rate of Reacovery')
plt.plot(dh_range, RoR_uneven_app * 100, linewidth = 2.5, color = 'gray', linestyle = '--', label = 'Rate of Reacovery (App.)')
plt.title('C. Rate of Recovery', loc = 'left', fontsize = 15)

plt.xticks(fontsize = 12.5); plt.xlabel('Step-up Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('%', fontsize = 15)
plt.legend(frameon = False)

plt.tight_layout()

#### Step down
dh_range_down = np.linspace(0, -0.05, 1000)
g = 9.81

# Compute v_post_opt
v_post_down_opt = np.sqrt(v_ave**2 + 2 * g * dh_range_down)

# Initialize a list to store the solutions
gama_opt_down_list = []

for v_post in v_post_down_opt:
    # Initial guess for gama (starting point for the solver)
    gama_guess = 0.35  # Adjust as necessary based on expected range for gama

    # Use fsolve to find the root of the equation
    gama_opt_down_solution = fsolve(equation_to_solve, gama_guess, args=(v_post, alpha_ave, v_ave))

    # Append the solution to the list
    gama_opt_down_list.append(gama_opt_down_solution[0])

# Convert the list to a NumPy array
gama_down_opt_array = np.array(gama_opt_down_list)

## Gama_opt must be less than 2 * alpha
gama_down_opt_array = gama_down_opt_array[gama_down_opt_array >= 0]

RoR_uneven_down     = 1 - (np.sin(2 * alpha_ave - gama_down_opt_array))**2 / (np.cos(gama_down_opt_array))**2
RoR_uneven_down_app = 1 - (2 * alpha_ave - gama_down_opt_array)**2

po_down = v_ave * np.tan(gama_down_opt_array); PO_down = 0.5 * po_down**2
co_down = v_ave * (np.sin(2 * alpha_ave - gama_down_opt_array)) / (np.cos(gama_down_opt_array)); CO_down = 0.5 * co_down**2

##
dh_range_down = dh_range_down[: len(gama_down_opt_array)]
x_range = [-0.04, -0.03, -0.02, -0.01, 0.00]
##
plt.figure(4, figsize=(14, 6)) 

plt.subplot(1,3,1)
plt.plot(dh_range_down, gama_down_opt_array, linewidth = 2.5, color = 'k')
plt.title('A. Mid-Transition Angle', loc = 'left', fontsize = 15)

plt.xticks(x_range, fontsize = 12.5); plt.xlabel('Step-down Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('Rad.', fontsize = 15)

plt.subplot(1,3,2)
plt.plot(dh_range_down, PO_down, label = 'Optimal Push-off', linewidth = 2.5, color = 'darkred')
plt.plot(dh_range_down, CO_down, label = 'Optimal Collision', linewidth = 2.5, color = 'navy')

plt.title('B. Step Transition Works', loc = 'left', fontsize = 15)
plt.legend(frameon = False)
plt.xticks(x_range, fontsize = 12.5); plt.xlabel('Step-down Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('J/kg', fontsize = 15)

plt.subplot(1,3,3)
plt.plot(dh_range_down, RoR_uneven_down * 100, linewidth = 2.5, color = 'k', label = 'Rate of Reacovery')
plt.plot(dh_range_down, RoR_uneven_down_app * 100, linewidth = 2.5, color = 'gray', linestyle = '--', label = 'Rate of Reacovery (App.)')
plt.title('C. Rate of Recovery', loc = 'left', fontsize = 15)

plt.xticks(x_range, fontsize = 12.5); plt.xlabel('Step-down Magnitude (m)', fontsize = 12.5)
plt.yticks(fontsize = 12.5); plt.ylabel('%', fontsize = 15)
plt.legend(frameon = False)

plt.tight_layout()

## average RoR for uneven walking
Uneven_RoR = np.append(RoR_uneven_down[::-1], RoR_uneven[: len(RoR_uneven_down)])
print('Average uneven walking RoR:', np.mean(Uneven_RoR))