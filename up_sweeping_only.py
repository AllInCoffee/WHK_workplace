# -*- coding: utf-8 -*-
"""
Created on Sun May  4 11:50:13 2025

@author: yanxi
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from numpy.linalg import eigvals
from scipy import constants
import math
import time

# Placeholder for physical constants (replace with actual values)
# feedback loop
feedback=0.2-0.25j


# constants under 1550nm
#Typical FCA Coefficient Value for Silicon (at 1550 nm): FCA cross section
sigma=1.45e-21 # m^2 
#The heat capacity per unit volume of silicon:
C= 1.63e6 # J/m^3
#Photon energy at 1550nm:
energy_ph=0.8*constants.elementary_charge #eV

# Kerr coefficient / SPM coefficient in terms of the nonlinear phase shift,
#gamma= 2*math.pi*n2/wavelength
gamma=3.1e-11 # m/W
# TPA-to-SPM strength ratio: r = beta/(2*gamma)
r=0.19
u = 30

# cavity constants
# Photon lifetime
tau_ph = 455e-12  # 455ps 
#Mode effective area 
A_eff=1e-13 # m^2 
T0=0 # initial temperature 

# free carrier lifetime depends on photon lifetime AND V_rb
# FC-lifetime decay constant
V_fc=2.15 # V
tau_c_sat=3.1 # 3.1ps! FC lifetime saturation
tau_c_0=140 # 140ps! FC lifetime when V_rb=0

# Thermal detuning coefficient 
delta_T= -9.7e9 # s^-1*K^-1
######################################################
# fixed TO decay time: 30ns
tau_th_absolute=30e-6 
tau_th_relative= tau_th_absolute/ tau_ph # near 66*tau_ph 

# group index
n_g=3.97
v_g= constants.c / n_g
#print(v_g)
# Ring-bus power coupling
theta=1.6e-3
# Round-trip time
t_R= 1.7e-12 # 1.7 ps
## input normalied factor
scale_in = t_R / math.sqrt(8 * tau_ph**3 * gamma * v_g * theta)
# print(scale_in) # 31998.6
## output normalied factor
scale_out = 1/math.sqrt(2*gamma*v_g*tau_ph)

scale_Delta= 1/ (2*tau_ph)

eta = 2 * delta_T * r / (gamma * v_g**2 * C)

#######################################################################
#################### Settings  ###########################################
Delta_MHz= -3000
tau_c_ps = 10
P_in=1e-10 # uW   ,  1e-10 ~ u_in_norminated 1e-6
#P_control_max= 25e-4 # W
P_control_max= 1e-3
num_pts =50


Delta = Delta_MHz * 1e6 *2 * math.pi/scale_Delta # norminated Delta
#u_in=math.sqrt(P_in*1e-6/A_eff) / scale_in


#u_control = math.sqrt(P_control/A_eff) / scale_in
if tau_c_ps > tau_c_0:
    # passive
    V_rb=0 
else:
    # active 
    V_rb = -V_fc * np.log((tau_c_ps - tau_c_sat) / (tau_c_0 - tau_c_sat))
tau_c_relative = tau_c_ps * 1e-12 / tau_ph
chi = (r * u * sigma) / (4 * energy_ph * gamma * v_g) * tau_c_relative
V_fact = 1 + V_rb * constants.elementary_charge / (2 * energy_ph)   

# Define your updated super_mode_real
def super_mode_real(t, S):
    x1, y1, x2, y2 = S
    s1 = x1 + 1j * y1
    s2 = x2 + 1j * y2

    s1_abs2 = x1**2 + y1**2
    s2_abs2 = x2**2 + y2**2

    ns1 = (chi * abs(s1+1j*s2)**4) / 4.0
    ns2 = (chi * abs(s2+1j*s1)**4) / 4.0
    Ts1 = eta * tau_th_relative * (
        V_fact * abs(s1+1j*s2)**4 / 4.0
        + ns1*abs(s1+1j*s2)**2 / (u * r * 2.0)
    ) + T0
    Ts2 = eta * tau_th_relative * (
        V_fact * abs(s2+1j*s1)**4 / 4.0
        + ns2*abs(s2+1j*s1)**2 / (u * r * 2.0)
    ) + T0

    ds1dt = -(1+1j*Delta) * s1 + 0.5j*(Ts1*(s1+1j*s2)+Ts2*(s1-1j*s2))\
        + 0.5*(1j-r)*(s1_abs2*s1+2*s2_abs2*s1-s2**2*(x1 - 1j * y1))\
        - 0.5*(1j+1/u)*(ns1*(s1+1j*s2)+ns2*(s1-1j*s2)) + u_in
    ds2dt = -(1+1j*Delta) * s2 + 0.5j*(Ts2*(1j*s1+s2)+Ts1*(s2-1j*s1))\
        + 0.5*(1j-r)*(s2_abs2*s2+2*s1_abs2*s2-s1**2*(x2 - 1j * y2))\
        - 0.5*(1j+1/u)*(ns2*(s2+1j*s1)+ns1*(s2-1j*s1))\
        + (feedback*s2 + u_control * 1j * np.sqrt(1 - abs(feedback)**2)) / (1 - feedback)

    return [ds1dt.real, ds1dt.imag, ds2dt.real, ds2dt.imag]

# Initial conditions to test
initial_conditions = [
    [0.0, 0.0, 0.0, 0.0],       # Base state (zero)
#     [0.1, 0.0, 0.0, 0.0],       # Perturbation in Re(s1)
#     [0.0, 0.1, 0.0, 0.0],       # Perturbation in Im(s1)
#     [0.0, 0.0, 0.1, 0.0],       # Perturbation in Re(s2)
#     [0.0, 0.0, 0.0, 0.1],       # Perturbation in Im(s2)
    
     [-0.1, 0.1, 0.1, 0.1],       # Equal phase excitation
#     # [-0.1, 0.1, -0.1, 0.1],     # Asymmetric phase flips
#     # [0.2, 0.0, 0.2, 0.0],       # In-phase excitation
#     # [0.0, 0.2, 0.0, -0.2],      # Quadrature phase shift
    
#     # [0.3, 0.3, 0.3, 0.3],       # Larger symmetric perturbation
#     # [-0.3, -0.3, 0.3, 0.3],     # Mixed-sign symmetry breaking
#     # [0.2, -0.2, -0.2, 0.2],     # Diagonal in state-space
#     # [0.4, 0.0, 0.0, 0.4],       # Excite Re(s1), Im(s2)
]
# #real_vals = np.linspace(-0.5, 0.5, 11)
# # #print(real_vals) # [-0.5 -0.4 -0.3 -0.2 -0.1  0.   0.1  0.2  0.3  0.4  0.5]
# real_vals = [-0.5,-0,3, -0.1, 0, 0.1, 0.3,0.5]
# # #imag_vals = np.linspace(-0.5, 0.5, 11)
# imag_vals = [-0.5,-0,3, -0.1, 0, 0.1, 0.3,0.5]
# initial_conditions = [(r, i, r, i) for r in real_vals for i in imag_vals]
# === Time Domain Parameters ===
t_span = [0, 20]
atol = 1e-9
rtol = 1e-8
# Jacobian approximation
def numerical_jacobian(func, x0, eps=1e-6):
    n = len(x0)
    jac = np.zeros((n, n))
    fx0 = np.array(func(x0))
    for i in range(n):
        dx = np.zeros(n)
        dx[i] = eps
        f1 = np.array(func(x0 + dx))
        jac[:, i] = (f1 - fx0) / eps
    return jac

# Refinement and stability checking
def refine_and_check(S_guess):
    sol = fsolve(lambda S: super_mode_real(0, S), S_guess)
    J = numerical_jacobian(lambda S: super_mode_real(0, S), sol)
    eigs = eigvals(J)
    is_stable = np.all(np.real(eigs) < 0)
    if not is_stable:
        print('not stable')
    return sol, is_stable

#Sweep function
# explore multistability (find all coexisting states at a given control power):
def bidirectional_sweep(P_list):
    P_vals, s1_abs2_list, s2_abs2_list, stabilities = [], [], [], []
    s1_list, s2_list=[],[]
    xi=[]
    for P_control in P_list:
        global u_control, u_in
        u_control = np.sqrt(P_control / A_eff) / scale_in
        u_in = np.sqrt(P_in / A_eff) / scale_in

        found = []
        for ic in initial_conditions:
            ans = solve_ivp(super_mode_real, t_span, ic, rtol=rtol, atol=atol)
            sol_guess = ans.y[:, -1]
            refined, is_stable = refine_and_check(sol_guess)

            if not any(np.linalg.norm(refined - np.array(prev)) < 1e-3 for prev in found):
                found.append(refined)
                s1 = refined[0] + 1j * refined[1]
                s2 = refined[2] + 1j * refined[3]
                s1_list.append(s1)
                s2_list.append(s2)
                P_vals.append(P_control)
                s1_abs2_list.append(abs(s1)**2)
                s2_abs2_list.append(abs(s2)**2)
                xi.append(np.sign(s1.real)*abs(s1))
                stabilities.append(is_stable)
    return P_vals, s1_list,s2_list,s1_abs2_list, s2_abs2_list, stabilities,xi
## note: but if we simulate hysteresis or dynamic evolution (e.g., forward/backward sweep): we need 
## to swap two "for" loops, it mimics a realistic experimental sweep, where the system "remembers" its previous state. It will help reveal hysteresis loops (i.e., path-dependence).
# Run bidirectional sweep
P_control_range = np.linspace(0, P_control_max, num_pts)
start=time.time()
P_up, s1_up, s2_up, s1_abs2, s2_abs2, stable_up,xi = bidirectional_sweep(P_control_range)

duration=time.time()-start
print(f'elapse: {duration/60:2f}min')

# Set global matplotlib parameters for consistency
# mpl.rcParams.update({
#     'font.size': 10,
#     'axes.labelsize': 10,
#     'axes.titlesize': 12,
#     'legend.fontsize': 10,
#     'xtick.labelsize': 6,
#     'ytick.labelsize': 6,
#     'figure.dpi': 300,
#     'savefig.dpi': 300,
#     'text.usetex': False  # Set to True if LaTeX is available
# })
############################################ xi and |s1|^2 stable branches
# Separate stable and unstable points
P_stable = [P for P, st in zip(P_up, stable_up) if st]
s1_stable = [s for s, st in zip(s1_up, stable_up) if st]
xi_stable = [xi for xi, st in zip(xi, stable_up) if st]
s1_abs2_stable=[s for s, st in zip(s1_abs2, stable_up) if st]
s2_abs2_stable=[s for s, st in zip(s2_abs2, stable_up) if st]

P_unstable = [P for P, st in zip(P_up, stable_up) if not st]
s1_unstable = [s for s, st in zip(s1_up, stable_up) if not st]
s1_abs2_unstable=[s for s, st in zip(s1_abs2, stable_up) if not st]
s2_abs2_unstable=[s for s, st in zip(s2_abs2, stable_up) if not st]

# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2)

# Plot stable branches for |s1|^2 on the first subplot (ax1)
ax1.plot(P_stable, s1_abs2_stable, 'o', color='blue', label=r'Stable $|s_1|^2$', markersize=4)
ax1.set_xlabel(r'$P_{control}$')
ax1.set_ylabel(r'$|s_1|^2$')
ax1.grid(True)
ax1.legend(loc='best')

# Plot stable branches for xi on the second subplot (ax2)
ax2.plot(P_stable, xi_stable, 'x', color='black', label=r'Stable $\xi$', markersize=4)
ax2.set_xlabel(r'$P_{control}$')
ax2.set_ylabel(r'$\xi$')
ax2.grid(True)
ax2.legend(loc='best')

ax1.tick_params(axis='both', which='major', labelsize=6)
ax2.tick_params(axis='both', which='major', labelsize=6)
ax1.legend(fontsize=10)
ax2.legend(fontsize=10)


fig.suptitle(r'Bifurcation Diagram: $|s_1|^2$ and $\xi$ vs Control Power', fontsize=14, y=0.93)
# Adjust layout to avoid overlap
plt.tight_layout()

plt.savefig("bifurcation_subplots_xi.png", dpi=300)  # Save to PNG file
# Show the plot
plt.show()

# Plot unstable branches
# ax.plot(P_unstable, s1_unstable, 'x', color='orange', label=r'Unstable $|s_1|^2$', markersize=4)
# ax.plot(P_unstable, s1_Re_unstable, 'x', color='red', label=r'Unstable Re{$s_1$}', markersize=4)
#################################### |s1|^2 and |s2|^2 ##########################
s1_abs2_stable=np.array(s1_abs2_stable)
s2_abs2_stable=np.array(s2_abs2_stable)
P_up = np.array(P_up)
xi = np.array(xi)
# Convert stability list to boolean array
stable_up = np.array(stable_up)

# Create masks
mask_stable = stable_up
mask_unstable = ~stable_up


# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2)

# --- Plot s1_up ---
ax1.scatter(P_up[mask_stable], s1_abs2_stable[mask_stable], 
            marker='o', color='blue', label='Stable $|s_1|^2$', s=10)
ax1.scatter(P_up[mask_unstable], s1_abs2_stable[mask_unstable], 
            marker='*', color='orange', label='Unstable $|s_1|^2$', s=10)
ax1.set_xlabel(r'$P_{control}$')
ax1.set_ylabel(r'$|s_1|^2$')
ax1.grid(True)
ax1.legend()

# --- Plot s2_up ---
ax2.scatter(P_up[mask_stable], s2_abs2_stable[mask_stable], 
            marker='o', color='green', label='Stable $|s_2|^2$', s=10)
ax2.scatter(P_up[mask_unstable], s2_abs2_stable[mask_unstable], 
            marker='*', color='red', label='Unstable $|s_2|^2$', s=10)
ax2.set_xlabel(r'$P_{control}$')
ax2.set_ylabel(r'$|s_2|^2$')
ax2.grid(True)
ax2.legend()

fig.suptitle(r'Bifurcation Diagram: $|s_1|^2$ and $|s_2|^2$ vs Control Power', fontsize=14, y=0.93)
plt.tight_layout()
plt.savefig("bifurcation_subplots_s1s2.png", dpi=300)  # Save to PNG file
plt.show()


import csv

# Save as CSV
with open('bifurcation_data.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['P_control (W)', 'Re(s1)', 'Im(s1)', '|s1|^2', 'Re(s2)', 'Im(s2)', '|s2|^2', 'Stable'])

    # Write each row
    for P, s1, s2, st, xi in zip(P_up, s1_up, s2_up, stable_up,xi):
        writer.writerow([
            P,
            s1.real, s1.imag, abs(s1)**2,
            s2.real, s2.imag, abs(s2)**2,
            int(st),  # 1 for stable, 0 for unstable
            xi
        ])
