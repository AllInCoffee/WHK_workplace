# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 07:34:54 2025

@author: yanxi
"""
from scipy.integrate import complex_ode
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
#import cmath
import math
import time
from sympy import *
from sympy import symbols
from scipy.optimize import fsolve

chi = 0.55
gamma_TPA=0.11
gamma_FCA=0.2
sigma_FCD=7.2
xi_T=0.074;
eta_lin=0.4;
eta_c=1
# set tau_th and tau_car
tau_car=10e-9
tau_th=100e-9
### set tau according to new tau_car
tau=1.85e9*tau_car
tao_theta=1.85e9*tau_th

## Delta= 100 (0,s_tilde_0) stable
### steady states:
### when s=0, s_tilde=s_tilde_0 , input for both MRR is half of P_control, |z|^2=|y|^2 s=0, super mode stable
Delta_fix=100
## [1.06,1.3] 
lower_limit = 1.001
upper_limit = 2.5
num_pts = 100
E_list = np.linspace(lower_limit, upper_limit, num_pts)
u_list=[]
for i in E_list:
    P_control = ((1 + gamma_TPA*i + gamma_FCA * tau*i**2)**2 + (Delta_fix + chi*i+ 
               tao_theta* xi_T*i*(eta_lin*eta_c+2*gamma_TPA*i + 2*gamma_FCA*tau*i**2) 
               -tau*i**2-sigma_FCD*tau**0.8*i**1.6)**2)*i
    u_list.append(math.sqrt(P_control/2))
plt.plot(u_list,E_list,label='|z|^2 (|y|^2)') 
# |z|^2=E=|s_tilde_0|^2/2

plt.title('Bifurcation Diagram of a Kerr Resonator with FCE and TOE')
plt.xlabel('$u_{control}$')
plt.ylabel('$|z_j|^2$')


## s_tilde_0 with feedback when s=0
plt.figure()
E=symbols('E')
Delta=symbols('Delta')
r=0.2-0.25j

# eqn2 = ((1- (r/(1-r)).real + gamma_TPA*E/2 + gamma_FCA * tau*E**2/4)**2 + (-(r/(1-r)).imag+Delta_fix + chi*E/2+ 
#            tao_theta* xi_T*gamma_TPA*E**2
#            -tau*E**2/4-sigma_FCD*(tau*E**2/4)**0.8)**2)*E*abs(1-r)**2/(1-abs(r)**2)
# u_list2=[]
# E_diff2=[]
# bi_u2=[]
# bi_E2=[]
# for j in E_list:
#     P_control=eqn2.subs([(E,j),(Delta,Delta_fix)])
#     u_list2.append(math.sqrt(P_control))
#     new=eqn2.diff(E).subs([(E,j),(Delta,Delta_fix)])
#     E_diff2.append(new)
# plt.plot(u_list2,E_list) ## u control vs s_tilde_0
# plt.title('Bifurcation Diagram of a Kerr Resonator with FCE and TOE')
# plt.xlabel('$u_{control}$')
# plt.ylabel('$|\\tilde{s}_0|^2$')

# ##### to find the working range , both s and s_tilde mode should be in non-bistable range so that it's in steady stable.
# zero_crossings = np.where(np.diff(np.sign(E_diff2)))[0]+1
# print(zero_crossings)
# if len(zero_crossings)==2:
#     bi_u2.append([u_list2[i] for i in zero_crossings])
#     bi_E2.append([E_list[i] for i in zero_crossings])
# print(bi_u2) #[[69.60141791858324, 9.281690767557794]]
# print(bi_E2) #[[1.060670707070707, 2.515167676767677]]
## so the range of E should be below 1, in this range any input u_-in would make the setting stable. 

###########################################################
## pertubate near steady state/ fixed point (0,s_tilde_0
## a small displacement from the fixpoint, assume s*==s, s_tilde constant


#### numerical 
# plt.figure()
#### 1. we can plot s but we can't extract xi from it, because psi and xi oscillate 
#### in original paper they assume ψ(t)relaxes much faster to its steady state than ξ(t)
def s_tilde_mode(t,s2,u_in,u_control):
    n0= tau*s2**4/4.0
    T0= tao_theta* xi_T*(s2)**2/2.0*(eta_lin*eta_c+gamma_TPA*(s2)**2 + gamma_FCA*2.0*n0)
    ### s_tilde mode steady state is when s == 0, so here ds2dt is not s1 dependent
    ds2dt= -(1- (r/(1-r)).real + gamma_TPA*abs(s2)**2/2.0 + gamma_FCA * n0)*s2-1j*s2*(Delta_fix + chi*abs(s2)**2/2.0+ 
                T0-n0-sigma_FCD*(n0)**0.8)+u_control*1j*math.sqrt(1-abs(r)**2)/(1-r)
    return ds2dt
def s_mode(t,s1,s2,n0,T0,u_in,u_control):
    ds1dt= -(1 + gamma_TPA*(abs(s1)**2+2*abs(s2)**2)/2.0 + gamma_FCA * n0)*s1-1j*s1* (Delta_fix + chi*(abs(s1)**2+2*abs(s2)**2)/2.0+ 
                T0-n0-sigma_FCD*(n0)**0.8  + 0.5*(1j*chi+gamma_TPA)*s2**2*np.conj(s1))+u_in
    return ds1dt

t0 = 0
tfinal = 100
ts = np.linspace(t0, tfinal, 200)

# u_in=0, u_control arbitray
u_in=0.1
u_control=2
parameter1=[0, u_control]
#### first, find the stable s_tilde when s = 0
s_tilde_stable=solve_ivp(s_tilde_mode,[t0, tfinal],[0],t_eval=ts,args=parameter1)
plt.plot(s_tilde_stable.t,s_tilde_stable.y[0],label='$\\tilde{S}_0$')
plt.title('$u_{{control}}={:.2f}$'.format(u_control),fontsize=15)
plt.legend()
plt.figure()
s2=s_tilde_stable.y[0,-1]  ## update s2 with stable value of s_tilde, it is real number
n0= tau*s2**4/4.0
T0= tao_theta* xi_T*(s2)**2/2.0*(eta_lin*eta_c+gamma_TPA*(s2)**2 + gamma_FCA*2.0*n0)
#print(s2,n0,T0)
parameter2=[s2,n0,T0,u_in, u_control]
s_initial=0+0j
ans=solve_ivp(s_mode,[t0, tfinal],[s_initial],t_eval=ts,args=parameter2)
plt.plot(ans.t, ans.y[0].real,label='$Re\{{S\}}$')
plt.plot(ans.t, ans.y[0].imag,label='$Im\{{S\}}$')
plt.xlabel('t')
plt.title('$s_{{initial}}={0}, u_{{in}}={1:.2f}, u_{{control}}={2:.2f}$'.format(s_initial,u_in,u_control),fontsize=15)
plt.legend()
psi=[math.atan(i/j)for i, j in zip(ans.y[0].imag,ans.y[0].real)]
plt.figure()
plt.plot(ans.t,[math.cos(k) for k in psi] )



###################################################################################

# plt.figure()
# for u_control in np.linspace(0.1, 5, 4):
#     for u_in in np.linspace(0.01, 0.1, 1):
#         parameter1=[u_in,u_control]
#         s_tilde_stable=solve_ivp(s_tilde_mode,[t0, tfinal],[0],t_eval=ts,args=parameter1)
#         s2=s_tilde_stable.y[0,-1]  ## update s2 with stable value of s_tilde, it is real number
#         n0= tau*s2**4/4.0
#         T0= tao_theta* xi_T*(s2)**2/2.0*(eta_lin*eta_c+gamma_TPA*(s2)**2 + gamma_FCA*2.0*n0)
#         ans=solve_ivp(s_mode,[t0, tfinal],[0],t_eval=ts,args=parameter2)
#         ##plt.plot(ans.t, ans.y[0].real,label='$Re\{{S\}}$, u_in={0:.3f},u_control={1:.2f}'.format(u_in,u_control))
#         ##plt.plot(ans.t, ans.y[0].imag,label='$Im\{{S\}}$, u_in={0:.3f},u_control={1:.2f}'.format(u_in,u_control))
#         plt.plot(ans.t, abs(ans.y[0]),label=' u_in={0:.3f},u_control={1:.2f}'.format(u_in,u_control))  ## abs(s) == |ξ|
#         ##plt.plot(ans.t, ans.y[1].real,label='$Re\{\tilde{S}\}$')
#         ##plt.plot(ans.t, ans.y[1].imag,label='$Im\{{\tilde{S}\}}$')
#         #plt.plot(ans.t, abs(ans.y[1])) ## oscilation without convergence

# plt.title('S Mode')
# #plt.axhline(y=0, color="grey", linestyle=":")
# plt.xlabel("t")
# plt.ylabel('$|\\xi|$')
# plt.ylim(0.07,0.09)
# plt.legend()

################################################################################
##### 2. transform s into a specific quadrature with xi and psi : s=xi*exp(psi*i)
#####    so that we can plot xi 

def spin(t,Y,s2,n0,T0,u_in,u_control):
    xi=Y[0]
    psi=Y[1]
    dpsidt=-(chi*(s2)**2+ (Delta_fix+xi_T*tau_th*gamma_TPA*s2**4)+chi*(xi)**2/2.0-0.5*math.cos(2*psi)*chi*s2**2\
                    +0.5*s2**2*gamma_TPA*math.sin(2*psi))+n0+sigma_FCD*n0**0.8
    dxidt=-(0.5*(2-chi*s2**2*math.sin(2*psi)+gamma_TPA*((xi)**2+2*(s2)**2-(s2)**2*math.cos(2*psi)))\
            +gamma_FCA*n0)
    
    return [dxidt,dpsidt]
print(s2,n0,T0)
plt.figure()
spin=solve_ivp(spin,[0, 100],[0,0],t_eval=np.linspace(t0, tfinal, 100),args=[s2,n0,T0,u_in,u_control])
print(spin)
plt.plot(spin.t, spin.y[0],label='$\\xi$')
#plt.plot(spin.t, spin.y[1],label='$\\psi$')
plt.legend()


# ####### s_tilde mode, always positive, no alternating signs
# plt.figure()
# for u_control in np.linspace(0.1, 5, 4):
#     for u_in in np.linspace(0.01, 1.1, 1):
#         parameter=[u_in,u_control]
#         ans=solve_ivp(super_mode,[t0, tfinal],y0,t_eval=ts,args=parameter)
#         plt.plot(ans.t, ans.y[1],label='$\\tilde{{s}}_{{mode}}$, u_in={0:.2f},u_control={1:.2f}'.format(u_in,u_control))

# plt.ylim(-0.005,0.1)
# plt.axhline(y=0, color="grey", linestyle=":")
# plt.xlabel("t")
# plt.title('$\\tilde{S}$ mode')
# plt.legend()
########### plot s vs u_control
# plt.figure()
# u_control_list=np.linspace(0.001, 10, 1000)
# result=[]
# for k in u_control_list: 
#     parameter=[0.05,k]
#     ans=solve_ivp(super_mode,[t0, tfinal],y0,t_eval=ts,args=parameter)
#     result.append(ans.y[0,-1]) # append the stabilized s
# plt.plot(u_control_list,result)
# plt.xlabel("$u_{control}$")
# plt.ylabel("S")




#####################
# ## found out the range of u for each ring to be bistable
# bi_u=[]
# bi_E=[]
# bi_delta=[]
# E=symbols('E')
# Delta=symbols('Delta') 
# eqn=((1 + gamma_TPA*E + gamma_FCA * tau*E**2)**2 + (Delta + chi*E+ 
#     tao_theta* xi_T*E*(eta_lin*eta_c+2*gamma_TPA*E + 2*gamma_FCA*tau*E**2) 
#     -tau*E**2-sigma_FCD*tau**0.8*E**1.6)**2)*E

# u_list=[]
# E_diff1=[]
# for j in E_list:
#     P=eqn.subs([(E,j),(Delta,Delta_fix)])
#     u_list.append(math.sqrt(P))
#     new=eqn.diff(E).subs([(E,j),(Delta,Delta_fix)])
#     E_diff1.append(new)
# zero_crossings = np.where(np.diff(np.sign(E_diff1)))[0]+1
# print(zero_crossings)
# if len(zero_crossings)==2:
#     bi_u.append([u_list[i] for i in zero_crossings])
#     bi_E.append([E_list[i] for i in zero_crossings])
# print(bi_u) #[[93.99630508619025, 9.52812341120595]]
# print(bi_E) #[[0.8182545454545455, 1.333388888888889]]
# #Delta=-100, there is bistability when z in in [93.99630508619025, 9.52812341120595]
# # |z|^2=E=|s_tilde_0|^2/2



