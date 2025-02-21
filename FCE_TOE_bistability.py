# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 22:27:24 2025

@author: yanxi
"""

from scipy.integrate import complex_ode
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
#import cmath
import math
import time

def CMT3(t,Y,Delta,u,chi,tau,k_T,gamma_TPA,gamma_FCA,sigma_FCD, xi_T,eta_lin,eta_c,tao_theta):
    z = Y[0]
    n = Y[1]
    T = Y[2]
    dzdt = -(k_T/2+1j*(Delta+chi*abs(z)**2))*z-math.sqrt(k)*u 
    +1j*(n+sigma_FCD*n**0.8-T)*z - (gamma_TPA*abs(z)**2+gamma_FCA*n)*z
    dndt = abs(z)**4-n/tau
    dTdt = xi_T*abs(z)**2*(eta_lin*eta_c + 2*gamma_TPA*abs(z)**2 + 2*gamma_FCA*n) - T/tao_theta
    return [dzdt, dndt, dTdt]

chi = 0.55
tau=18.5
k_l=50 * abs(chi)
k = 150 * abs(chi)
k_T=k+k_l
Delta_th= -math.sqrt(3)/2.0 * k_T
gamma_TPA=0.11
gamma_FCA=0.2
sigma_FCD=7.2

xi_T=0.074;
eta_lin=0.4;
eta_c=1;
tao_theta=185;

Delta=-100
u=np.arange(0,150,5)
y0 = [0+0j, 1+0j,1]
t_span=[0,1000]

forward=[]
backward=[]
start = time.time()
# from u=0 to high input power
for k in u: 
    parameter=[Delta,k,chi,tau,k_T,gamma_TPA,gamma_FCA,sigma_FCD, xi_T,eta_lin,eta_c,tao_theta]
    # get stablilized z and n
    ans=solve_ivp(CMT3,t_span,y0,args=parameter)
    forward.append(abs(ans.y[0,-1])**2) # append the stabilized |z|**2
end=ans.y[:,-1]
print(end)
# from high input power down to 0 
for k in list(reversed(u)): 
    parameter=[Delta,k,chi,tau,k_T,gamma_TPA,gamma_FCA,sigma_FCD, xi_T,eta_lin,eta_c,tao_theta]
    # get stablilized z and n
    ans=solve_ivp(CMT3,t_span,end,args=parameter)
    
    backward.append(abs(ans.y[0,-1])**2) # append the stabilized |z|**2

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Upward and Downward Integration')
ax1.plot(u, forward,'b',label='Upward: Delta={}'.format(Delta))  
ax2.plot(list(reversed(u)), backward,'r',label='Downward: Delta={}'.format(Delta))
ax1.legend(loc='upper left')
ax2.legend(loc='upper left')
fig.supxlabel('$u$')
fig.supylabel('$|z|^2$')

fig_new,ax3=plt.subplots()
fig_new.suptitle('Upward and Downward Integration')
ax3.plot(u, forward,'b',label='Upward: Delta={}'.format(Delta))  
ax3.plot(list(reversed(u)), backward,'r',label='Downward: Delta={}'.format(Delta))
ax3.legend(loc='upper left')
fig_new.supxlabel('u')
fig_new.supylabel('$|z|^2$')
end = time.time()
print('Time elapse: {:.1f} min'.format((end - start)/60))