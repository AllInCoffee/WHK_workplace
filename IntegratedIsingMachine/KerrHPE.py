# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:44:20 2025

@author: yanxi
"""

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
#import cmath
import math
import time
from sympy import *
from sympy import symbols
import numpy.ma as ma
from itertools import product

### global values
chi = 1## assume chi<0, delta >0, the opposite is also correct
k_l=50 * abs(chi)
k = 150 * abs(chi)
k_T=k+k_l

### analytical
lower_limit = 0.0001 # not zero to avoid nan
upper_limit = 400
num_pts = 200
E_list = np.linspace(lower_limit, upper_limit, num_pts)
E=symbols('E')
Delta=symbols('Delta')
eqn=((k_T/2)**2 + (Delta + chi*E)**2)*E/k
Delta_th= math.sqrt(3)/2.0 * k_T 
Delta_range=[i*Delta_th for i in [-0.7, -1.0, -1.3]]
#print (Delta_range)
P_diff=[]
############ plot bifurcation diagram illustrating bistability
### Method 1 --> overlapping problem, boarder points didn't disappear but still connected across missing points 
### solutions: https://matplotlib.org/stable/gallery/lines_bars_and_markers/masked_demo.html
# for j in Delta_range:
#     u_stable=[]
#     u_unstable=[]## temporary list to save all the u value for each Delta
#     E_stable=[]
#     E_unstable=[]
#     for i in E_list:
#         P=eqn.subs([(E,i),(Delta,j)])
#         temp=eqn.diff(E).subs([(E,i),(Delta,j)])
#         #print (temp, type(temp))
#         if temp <=0:
#             u_unstable.append(math.sqrt(P))
#             E_unstable.append(i)
#         else:
#             u_stable.append(math.sqrt(P))
#             E_stable.append(i)
#     plt.plot(u_stable,E_stable,'-',label='$\Delta = {:.1f}\Delta_{{th}} $'.format(j/Delta_th))
#     plt.plot(u_unstable,E_unstable,'--')
# plt.title('Bifurcation Diagram of a Single Kerr MRR')
# plt.xlabel('$u$')
# plt.ylabel('$|z|^2$')
# plt.legend()

# # ### Method 2 with mask--> color problem 
# # to imitate the color in original paper:
# my_color=['m','#1f77b4','#2ca02c']
# for ind,j in enumerate(Delta_range):
#     u_list=[]
#     u_stable=[]
#     u_unstable=[]## temporary list to save all the u value for each Delta
#     E_stable=[]
#     E_unstable=[]
#     # masks for points on statle and unstable branches
#     m=[]
#     for i in E_list:
#         P=eqn.subs([(E,i),(Delta,j)])
#         u_list.append(math.sqrt(P))
#         m.append( eqn.diff(E).subs([(E,i),(Delta,j)])<=0)
#     #print(m) # true values will be masked (not shown) ! 
#     mask_stable=m # the values of m must be copied over to new list, otherwise there is overlapping problem
#     mask_unstable=[not i for i in mask_stable]
#     plt.plot(np.ma.array(u_list,mask=mask_stable), np.ma.array(E_list,mask=mask_stable),'-',color=my_color[ind],label='$\Delta = {:.1f}\Delta_{{th}} $'.format(j/Delta_th))
#     plt.plot(np.ma.array(u_list,mask=mask_unstable), np.ma.array(E_list,mask=mask_unstable),'--',color=my_color[ind])
# plt.title('Bifurcation Diagram of a Single Kerr MRR')
# plt.xlabel('$u=\sqrt{P}$')
# plt.ylabel('$E=|z|^2$')
# plt.legend()



# ##########################
# ## choose a working point where |Delta|=1.3Delta_th
plt.figure()
Delta_op=-1.3*Delta_th
u_list_op=[]  ## if not renamed, remember to empty the u_list!
P_diff1=[]
for i in E_list:
    P=eqn.subs([(E,i),(Delta,Delta_op)])
    u_list_op.append(math.sqrt(P))
    # 1. order differential 
    new=eqn.diff(E).subs([(E,i),(Delta,Delta_op)])
    P_diff1.append(new)
# plot the curve E vs u under the working point Delta=-200 
# to see where there are bistable states
## plt.plot(u_list_op,E_list,label='$\Delta = {} $'.format(Delta_op))
#plt.plot(u_list_op,E_list,label='$\Delta =-225.17\chi $',color='#2ca02c')
####plot showing unstable branch
# mask_stable=[i<=0 for i in P_diff1] 
# mask_unstable=[not i for i in mask_stable]
# plt.plot(np.ma.array(u_list_op,mask=mask_stable), np.ma.array(E_list,mask=mask_stable),'-',color='#2ca02c',label='$\Delta =-225.17\chi $')
# plt.plot(np.ma.array(u_list_op,mask=mask_unstable), np.ma.array(E_list,mask=mask_unstable),'--',color='#2ca02c')

# plt.title('Bistability Range of u at a Fixed $\Delta = -1.3\Delta_{th}$')
# # plot the input range for bistability 
# # index where the Derivative change sign (approximately where f'(x)=0)
# zero_crossings = np.where(np.diff(np.sign(P_diff1)))[0]+1
# # find the boundries of u 
# u_range=[u_list_op[i] for i in zero_crossings]
# #print(u_range) # single resonator
# for i in u_range:
#     plt.axvline(x=i,color='red', linestyle='dashed')
    
# plt.xlabel('$u=\sqrt{P}$')
# plt.ylabel('$E=|z|^2$')
# plt.legend()




########### numerical part ######################

# def CMT1(t,z,Delta,u,chi,k_T,k):
#     dzdt = -(k_T/2+1j*(Delta+chi*abs(z)**2))*z-math.sqrt(k)*u
#     return dzdt

# Delta=-1.3*Delta_th
# u=200
# z0=[0+0j]
# t_span=(0,0.5)  ## from the plot we can see that it converges within 0.5

# ### plot time development of z, to see how long it takes to converge
# ### this time will be required in the t_span in numerical solver
# # parameter=[Delta,u,chi,tau,k_T,k]
# # ans=solve_ivp(CMT1,t_span,z0,args=parameter)
# # plt.plot(ans.t, ans.y[0].real,label='Re{z}')
# # plt.plot(ans.t, ans.y[0].imag,label='Im{z}')
# # plt.xlabel('t')
# # plt.legend()

# #### plot switching like bistability
# forward=[]
# backward=[]
# #start = time.time()
# u_list=np.linspace(0,301,300) # the more points the more steep 
# # from u=0 to high input power
# for i in u_list: 
#     parameter=[Delta,i,chi,k_T,k]
#     # get stablilized z and n
#     ans=solve_ivp(CMT1,t_span,z0,args=parameter)
#     #print('Upward stable z is: {}'.format(ans.y[0,-1]))
#     forward.append(abs(ans.y[0,-1])**2) # append the stabilized |z|**2
# end=ans.y[:,-1]
# #print(end)
# # from high input power down to 0, initial point is with high z and n, e.g. the end of upward or (100,100)
# for j in list(reversed(u_list)): 
#     parameter=[Delta,j,chi,k_T,k]
#     # get stablilized z and n
#     ans=solve_ivp(CMT1,t_span,end,args=parameter) 
#     #print(ans.y)
#     #print('Downward stable z is: {}'.format(ans.y[0,-1]))
#     backward.append(abs(ans.y[0,-1])**2) # append the stabilized |z|**2

# fig, (ax1, ax2) = plt.subplots(1, 2)
# fig.suptitle('Upward and Downward Integration')
# ax1.plot(u_list, forward,'b',label='Upward')  
# ax2.plot(list(reversed(u_list)), backward,'r',label='Downward')
# ax1.legend(loc="upper left")
# ax2.legend(loc="upper left")
# fig.supxlabel('$u=\sqrt{P}$')
# fig.supylabel('$E=|z|^2$') 
# plt.show()

# plt.figure()
# plt.title('Numerical Simulation of a Single Kerr MRR')
# plt.plot(u_list, forward,'b',label='Upward')  
# plt.plot(list(reversed(u_list)), backward,'r',label='Downward')
# plt.xlabel('$u=\sqrt{P}$')
# plt.ylabel('$E=|z|^2$')
# plt.legend()
# plt.show()

# finish = time.time()
# print('Time elapse: {:.1f} min'.format((finish - start)/60))
#################above is about single ring resonator #######################################
#################################################################################
########### super modes with feedback ############
chi = -1## assume chi<0, delta >0, the opposite is also correct
Delta=-225.17*chi
k_l=50 * abs(chi)
k = 150 * abs(chi)
k_T=k+k_l
r=0.2-0.25j
k_feedback=((1+r)/(1-r)).real*k
k_T_feedback= k_feedback+k_l
Delta_feedback=Delta+ (r/(1-r)).imag *k
print(k,k_T,k_feedback,k_T_feedback,Delta_feedback)
def s_tilde_mode(t,s2,u_in,u_control): ## s=0
    ds2dt= -(0.5*k_T_feedback+ 1j*Delta_feedback+1j*0.5*chi*s2**2)*s2 - math.sqrt(k_feedback)*1j*abs(1-r)/(1-r)*u_control
    return ds2dt
def s_mode(t,s1,s2,u_in): #
    ds1dt= -(0.5*k_T+ 1j*Delta+1j*0.5*chi*(abs(s1)**2+2*abs(s2)**2))*s1 + 0.5j*chi*s2**2*np.conj(s1)- math.sqrt(k)*u_in
    return ds1dt
#plt.figure()
# for s2 in range(0,10):
#     ans=solve_ivp(s_mode,[0,50],[0+0j],args=[s2,1]) ## !! s_initial must be complex, otherwise there is no change with change of s_tilde_0 
#     plt.plot(ans.t, ans.y[0])
def xi_ode(t,xi,s2):
    dxidt=-0.5*k_T*(1+(chi*s2**2/k_T)*math.sqrt(1-(2+2*Delta/(chi*s2**2)+(xi)**2/s2**2)**2))
    return dxidt

def super_mode(t,s,u_in,u_control):
    s1=s[0]
    s2=s[1]
    ds1dt= -(0.5*k_T+ 1j*Delta+1j*0.5*chi*(abs(s1)**2+2*abs(s2)**2))*s1 + 0.5j*chi*s2**2*np.conj(s1)- math.sqrt(k)*u_in
    ds2dt= -(0.5*k_T_feedback+ 1j*Delta_feedback+1j*0.5*chi*(abs(s2)**2+2*abs(s1)**2))*s2 + 0.5j*chi*s1**2*np.conj(s2) - math.sqrt(k_feedback)*1j*abs(1-r)/(1-r)*u_control
    return [ds1dt, ds2dt]

###### to get bistability through ode of s mode and simply mirror y axis #####################################################
# a=[]
# b=[]
# c=[]
# d=[]
# for u_control in np.linspace(100, 400,200):
#     parameters=[0.001,u_control]
#     ans=solve_ivp(super_mode, [0,50], [0+0j, 0+0j], args=parameters)
#     #print(ans.y[0,-1])
#     a.append(ans.y[0,-1].real)
#     b.append(ans.y[0,-1].imag)
#     c.append(abs(ans.y[0,-1]))
#     d.append(abs(ans.y[1,-1]))
# # print(ans.y[0,-1])
# # print(ans.y[1,-1])
# #plt.plot(np.linspace(100, 400,100),a,'g')
# #plt.plot(np.linspace(100, 400,100),b,'m')
# plt.plot([i/153.63 for i in np.linspace(100, 400,200)],c,'b')
# plt.plot([i/153.63 for i in np.linspace(100, 400,200)],[-i for i in c],'r')
# #plt.plot(np.linspace(100, 400,100),d,'r')

# low=-4*Delta/(3*chi)-math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# high=-4*Delta/(3*chi)+math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# print('y-axis low and high are: {0:.2f}, {1:.2f}'.format(low,high))
# ## y-axis low and high are: 204.31, 396.15
# for E2 in [low,high]:
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     #print(math.sqrt(P)) ## 143.98579683950538   177.7719851925709
#     plt.axvline(x=math.sqrt(P)/153.63,color='grey', linestyle='dashed')

# plt.ylabel('$\\xi$')
# plt.xlabel('$\\frac{u_{control}}{u_0}$')

#################### test if the relationship between u_control and s_tilde^2 is the same
####################      whether through super modes ode or calculating polynom
# a1=[]
# a2=[]
# b=[]
# ### input of ODE is u_control while input for polynom is E2 (s_tilde^2)
# ##ode
# for u_control in np.linspace(100, 400,10):
#     parameters=[0,u_control]
#     ans=solve_ivp(super_mode, [0,50], [0+0j,0], args=parameters) 
#     a1.append(abs(ans.y[1,-1])**2) ## s_tilde_0^2
    
#     parameters=[0.0001,u_control]
#     ans=solve_ivp(super_mode, [0,50], [0+0j,0], args=parameters) 
#     a2.append(abs(ans.y[1,-1])**2) ## s_tilde_0^2
# plt.plot(np.linspace(100, 400,10),a1, 'm',label='numerical:$u_{in}=0$')
# plt.plot(np.linspace(100, 400,10),a2, 'r',label='numerical:$u_{in}\\neq 0$')
# ## polynom
# ## the beginning and end of a1 and a2 is the same, so we can take a1 or a2 as reference
# ##for E2 in np.linspace(a2[0], a2[-1],10):  
# for E2 in np.linspace(a1[0], a1[-1],10):  
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     b.append(math.sqrt(P))
# plt.plot(b,np.linspace(a1[0], a1[-1],10),'g',label='analytical:$u_{in}=0$')
# ### two curves are different as long as u_in >0, egal wie klein, when u_in=0, they are the same
# ## as long as we give a bit u_in, there is bistability 
# plt.xlabel('$u_{control}$')
# plt.ylabel('$E_{cavity}$')
# plt.legend()


# ###### to get bistability through ODE of xi and psi Method 1 ##################
# plt.figure()
# x_axis=[]
# y_axis_plus=[]
# y_axis_minus=[]
# def xi_psi(t,Y,E2):
#     xi=Y[0]
#     psi=Y[1]
#     dxidt=-0.5*(k_T-chi*E2*math.sin(2*psi))*xi
#     dpsidt=-chi*E2-Delta-0.5*chi*xi**2+0.5*chi*E2*math.cos(2*psi)
#     return [dxidt, dpsidt]

# for E2 in np.linspace(100, 600,20):
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     #print(math.sqrt(P))
#     x_axis.append(math.sqrt(P))
#     #ans=solve_ivp(xi_ode,[0,50],[0],args=[E2]) ## math domain error
#     ans_plus=solve_ivp(xi_psi,[0,50],[0.001,0],args=[E2])
#     ans_minus=solve_ivp(xi_psi,[0,50],[-0.001,0],args=[E2])
#     y_axis_plus.append(ans_plus.y[0,-1])
#     y_axis_minus.append(ans_minus.y[0,-1])
# plt.plot(x_axis,y_axis_plus,'b')
# plt.plot(x_axis,y_axis_minus,'r')
# # for E2 in [low,high]:
# #     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
# #     plt.axvline(x=math.sqrt(P)/153.63,color='grey', linestyle='dashed')

# plt.ylabel('$\\xi$')
# plt.xlabel('$u_{control}$')

##################################################

# # ###### to get bistability through ODE of xi and psi Method 2: initial xi must be not zero ##################
########## more realistic than method 3 
plt.figure()
x_axis=[]
y_axis_plus=[]
y_axis_minus=[]
def xi_psi(t,Y,E2):
    xi=Y[0]
    psi=Y[1]
    dxidt=-0.5*(k_T-chi*E2*math.sin(2*psi))*xi
    dpsidt=-chi*E2-Delta-0.5*chi*xi**2+0.5*chi*E2*math.cos(2*psi)
    return [dxidt, dpsidt]
u_in=1e-10
for u_control in np.linspace(100, 200,100):
    #### practically, if u_in=0, xi can't be any number but zero. 
    #### we need initial non-zero xi to see bistability. so u_in can't be zero but must be small near to fixed point (0, s_tilde_0)
    #### there should be a relation between u_in and xi, because |xi|=|s|
    #### however, the change of u_in doesn't change much of the |xi|
    parameters=[u_in,u_control] 
    ans=solve_ivp(super_mode, [0,50], [0+0j,0], args=parameters)
    xi_abs=abs(ans.y[0,-1]) ## |s| == |xi|
    print(xi_abs)
    E2=abs(ans.y[1,-1])**2 ## s_tilde_0^2  ## diagram is in shape -->-- (dependent on supercritical or subcritical range??)
     
    ans_plus=solve_ivp(xi_psi,[0,50],[xi_abs,0],args=[E2])  ## condition of xi should be related to condition of u_in
    ans_minus=solve_ivp(xi_psi,[0,50],[-xi_abs,0],args=[E2])
    
    y_axis_plus.append(ans_plus.y[0,-1])
    y_axis_minus.append(ans_minus.y[0,-1])
plt.plot([i/153.63 for i in np.linspace(100, 200,100)],y_axis_plus,'b')
plt.plot([i/153.63 for i in np.linspace(100, 200,100)],y_axis_minus,'r')


low=-4*Delta/(3*chi)-math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
high=-4*Delta/(3*chi)+math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
print('y-axis low and high are: {0:.2f}, {1:.2f}'.format(low,high))
## y-axis low and high are: 204.31, 396.15
for E2 in [low,high]:
    P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
    #print(math.sqrt(P)) ## 143.98579683950538   177.7719851925709
    plt.axvline(x=math.sqrt(P)/153.63,color='grey', linestyle='dashed')

plt.ylabel('$\\xi$')
plt.xlabel('$\\frac{u_{control}}{u_0}$')

##################################################


# ###### to get bistability through ODE of xi and psi Method 3 -- final answer ##################
######### we get s_tilde_0 by u_in = 0, substitute s_tilde_0 into the ODE of psi and xi (these odes were given assumed u_in=0)
######### to see pertubation near (s,s_tilde)=(0,s_tilde_0), we have to pertubate xi a bit 
# plt.figure()
# x_axis=[]
# y_axis_plus=[]
# y_axis_minus=[]
# def xi_psi(t,Y,E2):
#     xi=Y[0]
#     psi=Y[1]
#     dxidt=-0.5*(k_T-chi*E2*math.sin(2*psi))*xi
#     dpsidt=-chi*E2-Delta-0.5*chi*xi**2+0.5*chi*E2*math.cos(2*psi)
#     return [dxidt, dpsidt]
# u_in=0
# for u_control in np.linspace(100, 200,500):
    
#     parameters=[u_in,u_control] 
#     ans=solve_ivp(super_mode, [0,50], [0+0j,0], args=parameters)
#     # xi_abs=abs(ans.y[0,-1]) ## |s| == |xi|
#     # print(xi_abs)
#     E2=abs(ans.y[1,-1])**2 ## s_tilde_0^2  ## diagram is in shape 
    
#     ### condition of xi should be related to condition of u_in
#     ## ans_plus=solve_ivp(xi_psi,[0,50],[xi_abs,0],args=[E2])  
#     ## ans_minus=solve_ivp(xi_psi,[0,50],[-xi_abs,0],args=[E2])
#     ans_plus=solve_ivp(xi_psi,[0,50],[0.0001,0],args=[E2]) 
#     ans_minus=solve_ivp(xi_psi,[0,50],[-0.0001,0],args=[E2])
#     y_axis_plus.append(ans_plus.y[0,-1])
#     y_axis_minus.append(ans_minus.y[0,-1])
# plt.plot([i/153.63 for i in np.linspace(100, 200,500)],y_axis_plus,'b')
# plt.plot([i/153.63 for i in np.linspace(100, 200,500)],y_axis_minus,'r')


# low=-4*Delta/(3*chi)-math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# high=-4*Delta/(3*chi)+math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# print('y-axis low and high are: {0:.2f}, {1:.2f}'.format(low,high))
# ## y-axis low and high are: 204.31, 396.15
# for E2 in [low,high]:
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     #print(math.sqrt(P)) ## 143.98579683950538   177.7719851925709
#     plt.axvline(x=math.sqrt(P)/153.63,color='grey', linestyle='dashed')

# plt.ylabel('$\\xi$')
# plt.xlabel('$\\frac{u_{control}}{u_0}$')



######################################### End ###########################################

### try to extract psi and xi out of s --> didn't work well
# psi=[math.atan(i/j)for i, j in zip(a,b)]
# cos_psi=[math.cos(i) for i in psi]
# x1=[i/j for i,j in zip(a,cos_psi)]
# sin_psi=[math.sin(i) for i in psi]
# x2=[i/j for i,j in zip(b,sin_psi)]
# plt.figure()
# #plt.plot([i/153.63 for i in np.linspace(100, 400,400)],x1,'r--',label='$\\xi$')
# #plt.plot([i/153.63 for i in np.linspace(100, 400,100)],x2,'b-',label='$\\xi$')
# sign_xi=np.zeros(len(a))
# for ind, value in enumerate(b):
#     if value==0:
#         #xi[ind]=a[ind]
#         sign_xi[ind]=np.sign(a[ind])
#     else:
#         tan_psi=a[ind]/b[ind]
#         if tan_psi>0:
#             sign_xi[ind]=np.sign(a[ind])
#         else:
#             sign_xi[ind]=-np.sign(a[ind])
            



#######################################################
# u_list=[]
# #P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*abs(s2)**2/2.0)**2)*s2**2/k_feedback
# for E2 in np.linspace(0,600,40):
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     u_list.append(math.sqrt(P))
# ### plot u_control vs s_tilde mode 
# plt.figure()
# plt.plot([i/153.63 for i in u_list],np.linspace(0, 600,40))
# plt.title('$\\tilde{S}$')
# low=-4*Delta/(3*chi)-math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# high=-4*Delta/(3*chi)+math.sqrt(4*Delta**2/(9*abs(chi)**2)-k_T**2/(3*abs(chi)**2))
# print('y-axis low and high are: {0:.2f}, {1:.2f}'.format(low,high))
# ## y-axis low and high are: 204.31, 396.15
# for E2 in [low,high]:
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     print(math.sqrt(P))
#     plt.axvline(x=math.sqrt(P)/153.63,color='red', linestyle='dashed')
###u_control for bistability is [143.98, 177.77]
# ####################################################
# plt.figure()
### note! different initial condition 0 or 0+0j have different results! 
# u_control=150
# parameter=[0,u_control]
# ans=solve_ivp(s_tilde_mode,[0,50],[0],args=parameter)
# s2=ans.y[0,-1]
# print(s2)
# ans=solve_ivp(s_tilde_mode,[0,50],[0+0j],args=parameter)
# s2=ans.y[0,-1]
# print(s2)
#######################################################
# parameter=[s2,1]
# ans=solve_ivp(s_mode,[0,50],[0+0j],args=parameter) ## initial point must be complex form! 
#plt.figure()
#plt.plot(ans.t, abs(ans.y[0]))
#print(u_list)
### plot u_control vs s mode 
# plt.figure()
# u_in=0.01
# xi0=[0]
# ## upward
# s_list=[]
# s_real=[]
# s_imag=[]
# xi_plus_list=[]
# xi_minus_list=[]
# u_range=[]
# for E2 in np.linspace(100, 396.15,100):
#     P=((k_T_feedback/2.0)**2 + (Delta_feedback+chi*E2/2.0)**2)*E2/k_feedback
#     u_range.append(math.sqrt(P))
# s2_low=math.sqrt(100)
# s2_high=math.sqrt(396.15)
# ## condition for dxidt = 0, xi**2=math.sqrt(s2**4-k_T**2/chi**2)-2*s2**2-2*Delta/chi
# #condition=math.sqrt(s2**4-k_T**2/chi**2)-2*s2**2-2*Delta/chi
# for s2 in np.linspace(s2_low,s2_high,100):
#     parameter=[s2,u_in]
#     ans=solve_ivp(s_mode,[0,50],[0+0j],args=parameter)
#     #print(ans.y[0,-1])
#     s_list.append(abs(ans.y[0,-1]))
#     s_real.append(ans.y[0,-1].real)
#     s_imag.append(ans.y[0,-1].imag)
#     condition1=s2**4-k_T**2/chi**2
#     #print(condition1)
#     #print(1-(2+2*Delta/(chi*s2**2))**2)
#     condition2=1-(2+2*Delta/(chi*s2**2))**2
#     if  condition2>0 and (condition1>0):
#         condition3=math.sqrt(condition1)-2*s2**2-2*Delta/chi
#         if condition3>0:
#             xi_plus_list.append(math.sqrt(condition3))
#             xi_minus_list.append(-math.sqrt(condition3))
#         else:
#             xi_plus_list.append(0)
#             xi_minus_list.append(0)
#     else:
#         xi_plus_list.append(0)
#         xi_minus_list.append(0)
#     # 
# ## plt.plot([i/153.63 for i in u_range], s_real, label='$Re\{S\}$')
# ## plt.plot([i/153.63 for i in u_range], s_imag, label='$Im\{S\}$')
# ## plt.plot([i/153.63 for i in u_range], s_list, label='$|S|$')
# #plt.plot([i/153.63 for i in u_range], xi_plus_list, label='$\\xi+$')
# #plt.plot([i/153.63 for i in u_range], xi_minus_list, label='$\\xi-$')

# plt.xlim(0.94,1.02)
# plt.ylabel('$\\xi$')
# plt.xlabel('$\\frac{u_{control}}{u_0}$')
# psi=[math.atan(i/j)for i, j in zip(s_imag,s_real)]
# cos_psi=[math.cos(i) for i in psi]
# x=[i/j for i,j in zip(s_real,cos_psi)]
# #plt.plot([i/153.63 for i in u_range],x,'r--',label='$\\xi$')
# plt.legend()





