# functions for Floquet analysis

#import h5py
import numpy as np
import math as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
#from scipy.stats import chi2
#from scipy import signal
#from scipy.fftpack import fft, fftshift
#import matplotlib.patches as mpatches
#from matplotlib.colors import colorConverter as cc
from datetime import datetime
import numpy.distutils.system_info as sysinfo
sysinfo.get_info('atlas')


# =============================================================================    
# time-step functions functions


def rk4_time_step( params , A , dt, stop_time , case_flag ):
  # uniform time step 4th-order Runge-Kutta time stepper

  time = 0. # non-dimensional time
  count = 0
  output_period = 10
  output_count = 0
  
  start_time_0 = datetime.now()
  while time < stop_time - dt: 

    #start_time_kcoeffs = datetime.now()
    k1 = rk4( params , time , A , count , 0 , case_flag )
    k2 = rk4( params , time + dt/2. , A + k1*dt/2. , count , 0 , case_flag )
    k3 = rk4( params , time + dt/2. , A + k2*dt/2. , count , 0 , case_flag )
    k4 = rk4( params , time + dt , A + k3*dt , count , 0 , case_flag )
    #time_elapsed = datetime.now() - start_time_kcoeffs
    #print('k coeff time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

    #start_time_Phi_update = datetime.now()
    A = A + ( k1 + k2*2. + k3*2. + k4 )*dt/6.; 
    #time_elapsed = datetime.now() - start_time_Phi_update
    #print('Phi update time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

    #output_count = perturb_monitor( time , count , output_count , output_period , Phin , z , params , 'plot' )
    time = time + dt # non-dimensional time
    count = count + 1
    #freq = 100
    if params['freq'] != 0:
        if np.floor(count/params['freq']) == count/params['freq']:
            print( '%.2f complete' %(count/params['Nt']) )
            time_elapsed = datetime.now() - start_time_0
            print('Wall time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
    check_matrix(Phin,'Phin')

  dtf = stop_time - time
  k1 = rk4( params , time , A , count , 0 , case_flag )
  k2 = rk4( params , time + dtf/2. , A + k1*dtf/2. , count , 0 , case_flag )
  k3 = rk4( params , time + dtf/2. , A + k2*dtf/2. , count , 0 , case_flag )
  k4 = rk4( params , time + dtf , A + k3*dtf , count , 0 , case_flag )
  A = A + ( k1 + k2*2. + k3*2. + k4 )*dtf/6.; 

  # this is where conservation of mass needs to be checked

  final_time = time + dtf # non-dimensional final time
  count = count + 1

  if params['freq'] != 0:
      print( '%.2f complete' %(count/params['Nt']) )
      print('Wall time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

  #print('RK4 method, final time = ', final_time)
  return A, final_time


def rk4( params , time , A , count , plot_flag , case_flag ):
  # 4th-order Runge-Kutta functions 
  
  if case_flag == 'toy':
      for j in range(0,params['Ng']):
          if j == params['m']:
              krk = -params['nu']*n[j]**2.*A[j] + 1.
          else:
              krk = -params['nu']*n[j]**2.*A[j]

  return krk

