#
# Bryan Kaiser
# /2019



#import h5py
import numpy as np
import math as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import scipy
#from scipy import signal
#import sys
#sys.path.insert(0, '/path/to/application/app/folder')
#import functions as fn

import functions as fn
from datetime import datetime

figure_path = "./figures/"


# =============================================================================


Nz = 100 # number of points in z
z = np.linspace(0.,np.pi,num=Nz,endpoint=True)
H = np.pi
u0 = z*(np.pi*np.ones(np.shape(z))-z)

Nf = 100 # number of Fourier modes
n = np.linspace(1.,Nf,num=Nf,endpoint=True)
u0f = np.zeros([Nz])
for j in range(0,Nz):
    u0f[j] = sum(4./np.pi*(1.-(-1.)**n)/(n**3.)*np.sin(n*z[j]))



Nt = 1000
t = np.linspace(0.,2.*np.pi,num=Nt,endpoint=True)
m = 5.
nu = 1e-3

u = np.zeros([Nz,Nt])
for i in range(0,Nt):
    for j in range(0,Nz):
        u[j,i] = (1.-np.exp(-nu*m**2.*t[i]))*np.sin(m*z[j])/(nu*m**2.) + \
        sum(4./np.pi*(1.-(-1.)**n)/(n**3.)*np.sin(n*z[j])*np.exp(-nu*n**2.*t[i]))


T,Z = np.meshgrid(t/(2.*np.pi),z/H)
#print(np.shape(Z))

plotname = figure_path +'analytical_solution.png' 
plottitle = r"$u(z,t)$ " 
fig = plt.figure(figsize=(16, 8))
CS = plt.contourf(T,Z,u,200,cmap='seismic')
plt.xlabel(r"t",fontsize=16);
plt.ylabel(r"z",fontsize=16); 
plt.colorbar(CS)
plt.title(plottitle,fontsize=16);
plt.savefig(plotname,format="png"); plt.close(fig);




# Galerkin method

dt = t[1]-t[0]
#print(dt)
#print(t[0])
Ng = 7 # number of Galerkin modes
n = np.linspace(1.,Ng,num=Ng,endpoint=True)

# A matrix & initial condition
A = np.zeros([Ng,Nt])
for j in range(0,Ng):
    A[j,0] = 4.*(1.-np.cos(n[j]*np.pi))/(np.pi*n[j]**3.)
   
#print(A[10,0])

def rk4( params , A ):
    #print(A)
    #print(n)
    krk = np.zeros([Ng])
    # 4th-order Runge-Kutta functions  
    for j in range(0,params['Ng']):
          if j == (params['m']-1):
              krk[j] = -params['nu']*n[j]**2.*A[j] + 1.
          else:
              krk[j] = -params['nu']*n[j]**2.*A[j]
          #if j == 10:
          #    print(-params['nu']*n[j]**2.*A[j])
    #print(krk)
    return krk

params = {'Ng':Ng, 'nu':nu, 'm':m}

for k in range(1,Nt):
    #print(t[k-1])
    k1 = rk4( params , A[:,k-1] )
    k2 = rk4( params , A[:,k-1] + k1*dt/2. )
    k3 = rk4( params , A[:,k-1] + k2*dt/2. )
    k4 = rk4( params , A[:,k-1] + k3*dt )
    A[:,k] = A[:,k-1] + ( k1 + k2*2. + k3*2. + k4 )*dt/6.



#print(np.shape(n))
ug = np.zeros([Nz,Nt])
for k in range(0,Nt):
    for j in range(0,Nz):
        ug[j,k] = sum(np.multiply(A[:,k],np.sin(n[:]*z[j])))


plotname = figure_path +'Galerkin_solution.png' 
plottitle = r"$u(z,t)$ " 
fig = plt.figure(figsize=(16, 8))
CS = plt.contourf(T,Z,ug,200,cmap='seismic')
plt.xlabel(r"t",fontsize=16);
plt.ylabel(r"z",fontsize=16); 
plt.colorbar(CS)
plt.title(plottitle,fontsize=16);
plt.savefig(plotname,format="png"); plt.close(fig);


plotname = figure_path +'initial_condition.png' #%(start_time,end_time)
plottitle = r"initial condition" #, $\tau_w/U_0^2$ Re=%.2f, Pr=%.1f" #%(Re,Pr)
fig = plt.figure(figsize=(6,6))
plt.plot(u0,z,'b',label=r"$z(\pi-z)$")
plt.plot(u0f,z,'--r',label=r"Fourier")
plt.plot(u[:,0],z,'--g',label=r"Fourier")
plt.plot(ug[:,0],z,'or',label=r"Galerkin")
plt.xlabel(r"$u_0$",fontsize=13); 
plt.ylabel(r"z/H",fontsize=13); 
plt.title(plottitle);
plt.legend(loc=1,fontsize=13)
plt.grid()
plt.savefig(plotname,format="png"); plt.close(fig);
