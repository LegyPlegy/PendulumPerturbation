# -*- coding: utf-8 -*-
"""
Pendulum_Perturbation.py is a script that does analysis on a double pendulum,
specifically on the divergence in trajectories for double pendulums with similar
but slightly varying initial conditions.

Team Lead: Jorge Ramirez
Engine code: Angel Torres
Analysis: Maia Adams

Supervised by Dan Lathrop as a part of the University of Maryland, College Park
PHYS410: Classical Mechanics

"""

"""
This code is split into three sections.

    [0] simulation settings
        a chunk of code dedicated to simulation settings that are global
        throughout the entire script
        
    [1] numerical computation
        a chunk of code dedicated to the numerical computation of 
        double pendulum trajectories    
    
    [2] analytics
        a chunk of code dedicated to running analytics on the computations,
        e.g. lyapunov exponent, chaos measurements, etc
        
    [3] plotting algorithm
        a chunk of code dedicated to creating visual representations of the
        various double pendulum trajectories in an aesthetically-pleasing
        realtime plot format.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sp
from time import time



'''
===============================================================================
[0] simulation settings
    Inputs -
        user defined constants
        universal constants
    Outputs - 
        global variables to be used in [1] [2] (3)
===============================================================================        
'''

# first declare our simulation variables

m1, m2 = 1, 1  # pendulum mass 1 and mass 2  [kg]
l1, l2 = 1, 1  # pendulum length 1 and length 2 [meters]
g = 9.81  # gravitational acceleration

global_constants = (g, m1, m2, l1, l2)




'''
===============================================================================
[1] numerical computation (Angel Torres)
    Inputs - 
        time series 
        global settings
        
    Outputs -
        array of computed trajectory from section
===============================================================================     
'''

def double_pen(z, t, param):
    """
    Defines the diff eq for the double pendulum
    
    Inputs:
        z - list of state variables
            [theta1_dot, theta1_dotdot, theta2_dot, theta2_dotdot]
            
        t - time series
        
        param - global variables
            [m1, m2, l1, l2, g]
        
    Outputs:
        f - function f(t1dot, t1dotdot, t2dot, t2dotdot)
                
    """
    theta1, T1, theta2, T2 = z
    m1, m2, l1, l2, g = param
    C = np.cos(theta1 - theta2)
    c2 = np.cos(2*theta1 - 2*theta2)
    S = np.sin(theta1-theta2)
    M = m1+m2
    
#Simplifications for repetition between ODES for theta1doubledot and theta2doubledot 
#Four 1st Order Equations for the Coupled Second Order ODEs
    
    T1d = (-g*(2*m1+m2)*np.sin(theta1)-m2*g*S-2*S*m2*((T2**2)*l2+T1**2*l1*C))/(l1*(2*m1+m2-m2*c2))
    T2d = (2*S*((T1**2)*l1*M+g*M*np.cos(theta1)+(T2**2)*l2*m2*C))/(l2*(2*m1+m2-m2*c2))
    
    f = [T1, T1d, T2, T2d]
    
    return f

#Global Parameters
#Mass
m1 = 1
m2 = 1

#Lengths
l1 = 1
l2 = 1

#gravity
g = 9.8
def simulate_pendulum(time_series, init_cond, global_settings):
    """
  
    a method that simulates the trajectory of ONE pendulum for a given
    init cond and time series
        
        Inputs - 
            time_series: 
                [float array] a numpy array containing time steps
                created using numpy.linspace
              
            init_cond: 
                [tuple] a tuple of the form (Є1, Є2) where 
                Є1 is shift in θ1 and Є2 is shift in θ2
                
            global_settings: 
                [tuple] a tuple containing the simulation parameters
                        see section [0] for more info
  
        Outputs -
            sim_trajectory:
                [2xn array] a 2D array where column 1 and 2 contain the final 
                trajectory for θ1 and θ2, where n is the total steps
                
    """
    
        
    sim_trajectory = [0, 1, 2]  # contains computed values for f(θ1) and f(θ2) 
    
    return sim_trajectory
            




'''
===============================================================================
[2] analytics (Maia Adams)
    Inputs -
        time series
        array of computed trajectory from section [1]
        
    Outputs - 
        lyapunov exponent
        evolution of chaos??
        divergence of phase-space curves??
===============================================================================
'''

def lyapunov_exp(sim_trajectory, time):
    
    delta = np.abs(sim_trajectory[1]-sim_trajectory[0]) 
    # takes the difference between trajectories f(θ1) and f(θ2)
    
    lyapunov = np.log(delta)/time #lyapunov exponent
   
    return lyapunov




'''
===============================================================================
(3) plotting algorithms (Jorge Ramirez)
    Inputs - 
        time series
        array of computed trajector from sectino [1]
    
    Outputs -
        nxn plot of realtime pendulums
===============================================================================
'''

def init():
    """initialize animation"""
    line.set_data([], [])  # start with empty sets
    time_text.set_text('')  # make sure time is 0 at the beginning
    return line, time_text

def animate(i):
    """used for FuncAnimation: this is iterated over i frames"""
    
    # define function and store in xdata and ydata, iterating over i
    t = 0.1*i
    x = t*np.sin(t) 
    y = t*np.cos(t) 
    xdata.append(x) 
    ydata.append(y) 
    
    # update the line with new data
    line.set_data(xdata, ydata) 
    
    # keep track of time
    time_text.set_text("time = {0:.2f}s".format(t) )
    
    return line, time_text

# create figure 
fig = plt.figure() 
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-50, 50), ylim=(-50, 50))  
ax.grid()  # add grid to figure

# now, begin our animation
line, = ax.plot([], [], lw=2)  # required for FuncAnimation
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes) # add text to top left

# begin the animation method
xdata, ydata = [], []  # store data in here
anim = animation.FuncAnimation(fig, animate, init_func=init, 
							frames=500, interval=20, blit=True) 

# save animation
anim.save('double_pendulum.gif') 
print("Done")


