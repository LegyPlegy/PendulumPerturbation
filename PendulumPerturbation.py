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


# amount to dither the initial conditions by

perturbation = 0.005

#num_frames = 100
num_pendulums = 7

# Initial Conditions
# theta1 and theta2 assign the inital angle of the double pendulum, degrees
theta1 = 70 - (num_pendulums*perturbation*0.5)  # have 90 as the center 
theta2 = 0
#T1 and T2 are the initial velocities
T1 = 0.0
T2 = 0.0

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

t = np.linspace(0, 25, 500)
y0 = [theta1, theta2, T1, T2]
z = sp.integrate.odeint(double_pen, y0, t, args = (global_constants))




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

def lyapunov_exp(ref_pend, crazy_pend, num_pend, time):
    '''
    Inputs:    
        ref_pend - trajectory of pendulum with baseline init 
        crazy_pend - pendulum with highest perturbation (ref_pend + n*perturbation)
        num_pend - the number of pendulums
        time - time series for plotting
        
    Outputs:
        lyapunov exponent
    '''

    ref_pend_dist = np.sqrt(ref_pend[0]**2 + ref_pend[1]**2)
    crazy_pend_dist = np.sqrt(crazy_pend[0]**2 + crazy_pend[1]**2)
    
    delta = np.abs(ref_pend_dist - crazy_pend_dist) 
    
    if time == 0: # do not divide by zero!
        time += 0.00001
    
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

# create figure 
fig = plt.figure() 
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
                     xlim=(-2, 2), ylim=(-2, 2))  
ax.grid()  # add grid to figure
plt.title("Lower Half of Double Pendulum: Trajectory Analysis")


# first, create a series of line objects and store them
lines = []
for index in range(num_pendulums):
    lobj = ax.plot([],[],lw=2, alpha=0.5)[0]
    lines.append(lobj)

# call our method to make some pendulum trajectories
trajectories = []
for num in range(num_pendulums):
    y0[0] += perturbation  # increase theta 1 by 0.1 rad every time
    
    # generate trajectories using odeint
    z = sp.integrate.odeint(double_pen, y0, t, args = (global_constants))
    
    # unpack our theta 1 and theta 2
    theta1, theta2 = z[:,0], z[:,2]
    x1 = global_constants[3] * np.sin(theta1)
    y1 = -global_constants[4] * np.cos(theta1)
    
    x2 = x1 + global_constants[4] * np.sin(theta2)
    y2 = y1 - global_constants[4] * np.cos(theta2)
    trajectories.append((x2, y2))
    
# create text objects for 
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes) # add text to top left
lyapunov_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)

# create an initializer to empty every line object
def init():
    for line in lines:
        line.set_data([],[])
    return lines

# create a list of lists to keep track of data
all_xsets, all_ysets = [], []
for line in lines:
    all_xsets.append([])
    all_ysets.append([])

# create our animation function, that updates every line object every frame
def animate(i, data_list, all_xsets, all_ysets, lines, time_text, lyapunov_text):
    """
    Inputs:
        i - variable used to animate in FuncAnimate
        data_list - list of tuples with each tuple having a x and y array 
        lines - list of line objects
        
    Outputs:
        idk, it's used inside FuncAnimate
    """
    
    # update the dataset for each line object dynamically
    for line_num, line in enumerate(lines):
        # extract data from data-list
        next_xval = data_list[line_num][0][i]
        next_yval = data_list[line_num][1][i]
        
        # store that into each line's designated "all_xset" and "all_yset"
        all_xsets[line_num].append(next_xval)
        all_ysets[line_num].append(next_yval)
        
        # get rid of old traces to declutter graph
        if len(all_xsets[line_num]) >= 200:
            del all_xsets[line_num][0:50]
            del all_ysets[line_num][0:50]
        
                
        # finally, dynamically update that line's data with extracted values
        line.set_data(all_xsets[line_num], all_ysets[line_num])
        
        # calculate lyapunov exponent every frame
        # use furthest pendulum and the base pendulum
        
        if line_num == (num_pendulums-1)/2:  # base pendulum
            base_pendulum = (next_xval, next_yval)
        if line_num == num_pendulums-1:
            crazy_pendulum = (next_xval, next_yval)
        
    lyapunov = lyapunov_exp(base_pendulum, crazy_pendulum, num_pendulums, t[i])
    
    # update the text every frame
    time_text.set_text("Time [s]: {:.2f}".format(t[i]))
    lyapunov_text.set_text("Largest Lyapunov: {:.3f}".format(lyapunov))
    
        
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, interval=38,
                               save_count=sys.maxsize,
                               fargs=(trajectories, all_xsets, all_ysets, 
                                      lines, time_text, lyapunov_text))




# save animation
anim.save('double_pendulum.gif') 
plt.show()



