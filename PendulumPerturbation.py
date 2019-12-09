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
from scipy import integrate


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


num_frames = 100


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

def double_pen(z, t, m1, m2, l1, l2, g):
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

# Initial Conditions
# theta1 and theta2 assign the inital angle of the double pendulum
theta1 = 0.1
theta2 = 0.1
#T1 and T2 are the initial velocities
T1 = 0.0
T2 = 0.0
t = np.linspace(0, 50, 501)
y0 = theta1, theta2, T1, T2
z = sp.integrate.odeint(double_pen, y0, t, args = (global_constants))
'''
 
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
                
'''             




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


def animate(i, data_list, lines):
    """used for FuncAnimation: this is iterated over i frames"""
    
    
    for num, line in enumerate(lines):        
        # update the line with new data
        lines[num][0].set_xdata(data_list[num][0][:i])
        lines[num][0].set_ydata(data_list[num][1][:i])
        
    # keep track of time
    #time_text.set_text("time = {0:.2f}s".format(t))
    
    return lines[:]


# create figure 
fig = plt.figure() 
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-10, 10), ylim=(-10, 10))  
ax.grid()  # add grid to figure

# create fake data
omegas, trajectories = [2, 5, 9, 6], []
t = np.linspace(0, 5*np.pi, num_frames)

for omega in omegas:
    # trajectories is a list of 2 arrays, the first being sinAt and second being cosAt
    trajectories.append((t*np.sin(omega*t), t*np.cos(omega*t)))

# create a set of line objects for all our datasets

trajectory_lines = []
for trajectory in trajectories:
    trajectory_lines.append(ax.plot([], [], lw=1))  # required for FuncAnimation


#time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes) # add text to top left

# init the animation method 
anim = animation.FuncAnimation(fig, animate, fargs=[trajectories, trajectory_lines], \
                               frames=num_frames, interval=20, blit=True) 

# save animation
#anim.save('double_pendulum.gif') 
print("Done")



