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



'''
[0] simulation settings
    Inputs -
        user defined constants
        universal constants
    Outputs - 
        global variables to be used in [1] [2] (3)
        
'''



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
import numpy as np
import matplotlib as plt
import scipy as py

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
t = np.linspace(0, 5, 0.1)
traj = []

def lyapunov_exponent()



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
