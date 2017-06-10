# Galerkin spectral method p.47 Zikanov (2010)
# Bryan Kaiser
# 6/19/16

import math as ma
import numpy as np
import pylab as pl
import cmath as cma

# Numerical Galerkin spectral solution to:
# du/dt = a^2*(d^2u/dx^2)+sin(5x) 
# in the domain 0 < x < pi for 0 < t < T and the boundary conditions:
# u(0,t)=u(pi,t)=0 
# and the initial condition:
# u(x,0)=x(pi-x)

# note: use ''' and ''' again for multiline comments.


#==============================================================================
# functions

'''
def RK4(An,nu,t):
for n in range(1,Nc):
    if n != 5:
        for i in range(0,Nx-1):
        	An[n][i] = ((4.0*(1.0-ma.cos(np.pi*n)))/(np.pi*n**3.0))
	k
	return k
'''


#==============================================================================
# set up parameters

# domain / time series
Nc = 100 # number of coeffs
Nx = 60 # number of spatial points
Nt = 1000 # timesteps
T = 3000.0 # s, final time
t = 4.0 # s, initial time
dt = T/Nt # s, time step

# physical
a = np.sqrt(0.5) # sqrt(m^2/s), diffusion constant
x = np.linspace(0,np.pi,Nx) # m
A0 = np.zeros((Nc,Nx)) # An(t)*sin(n*x) # An(t)*sin(n*x), 
An = np.zeros((Nc,Nx)) # An(t)*sin(n*x) # An(t)*sin(n*x), 
A5 = np.zeros((1,Nx)) # A5(t)*sin(5*x) # A5(t)*sin(5*x), forcing term
uic = [x[i]*(np.pi-x[i]) for i in range(len(x))] # m/s, initial u(x,0)

#print a
#print np.shape(An)
#print np.shape(A5)

#==============================================================================
# compute initial coefficients An for An(t)*sin(n*x) at t = 0:

for n in range(1,Nc):
    if n != 5:
        for i in range(0,Nx-1):
        	A0[n][i] = ((4.0*(1.0-ma.cos(np.pi*n)))/(np.pi*n**3.0))

# now RK4 advance An from A0




#==============================================================================
# compute coefficients An for An(t)*sin(n*x):

t=T

for n in range(1,Nc):
    if n != 5:
        for i in range(0,Nx-1):
        	An[n][i] = ((4.0*(1.0-ma.cos(np.pi*n)))/(np.pi*n**3.0))*\
		ma.exp(-t*(n*a)**2.0)*ma.sin(n*x[i])
    elif n == 5:
        for i in range(0,Nx-1):
        	A5[0][i] = (((4.0*(1.0-ma.cos(np.pi*5.0)))/(np.pi*5**3)\
		-np.pi/(50*a**2))*ma.exp(-t*(5.0*a)**2.0)+ma.pi/(50.0*a**2))\
		*ma.sin(5.0*x[i]);


#==============================================================================
# solution

u = sum(An)+A5;
u = np.reshape(u, (1,np.product(u.shape)))[0]

pl.plot(x,uic,'k',x,u,'b')
pl.ylabel('u')
pl.xlabel('x')
#pl.xscale('log')
pl.show()










#==============================================================================
# DECLARED FUNCTIONS

def dft(f,dt):
	N = len(f) # length of signal
	ii = [ii**1 for ii in range(0,N)] # Hz
	Fs = 1.0/dt # s^{-1} 
	F = np.zeros(shape=(N,N),dtype=complex)
	for n in range(0,N-1):
		for k in range(0,N-1):
			F[n][k]=complex(f[k]*cma.exp(-1j*(2*ma.pi*n/N)*k))
	F = np.sum(F,axis=1,dtype=complex)
	for i in range(0,N-1):
		F[i] = abs(F[i])*(2.0/N)
	F = F.real
	#F[0] = F[0]/2.0 # dc correction, output
	hz = []
	for k in range(0,N):
        	hz.append(ii[k]*Fs/N)
	return F, hz

#==============================================================================

#a = 1j
#print a**2

# The sampling rate and time domain parameters:
Fs = 128 # samples/s, sampling rate
dt = 1.0/Fs # s, sampling times
N =2*Fs # number of samples, must be at least 2*Fs
k = [k**1 for k in range(0,N)] # sample index
t = []
for i in range(0,len(k)):
	t.append(k[i]*dt) # s, discrete sampling times
T = t[N-1] # s, total sequence time

# The sampled signal:
f = np.zeros(N) # signal f(t)
for i in range(N): # (0 to N-1)
	f[i] = 5+ma.cos(t[i]*ma.pi*2.0)*3.0+ \
	ma.cos(t[i]*ma.pi*4.0)*5.0+ \
	ma.cos(t[i]*ma.pi*6.0)*7.0

# Sampled signal plot:
#pylab.ylabel('f(t)')
#pylab.xlabel('t')
#pylab.show(pylab.plot(t,f))

#==============================================================================

# python built-in Fast Fourier Transform (FFT)
F = np.fft.fft(f) # unscaled frequency domain F(omega)
N = 2*Fs*1.0 # switch to floating point value
for i in range(0,len(F)):
	F[i] = abs(F[i])*(2/N)
hz = []
for i in range(0,len(k)):
        hz.append(k[i]*Fs/N)
pl.ylabel('F')
pl.xlabel('hz')
pl.xscale('log')
pl.show(pl.plot(hz,F))

# my dft function:
[F2,hz2] = dft(f,dt)
pl.ylabel('F2')
pl.xlabel('hz')
pl.xscale('log')
pl.show(pl.plot(hz2,F2))


