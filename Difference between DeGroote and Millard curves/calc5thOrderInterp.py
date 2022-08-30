# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calc5thOrderInterp.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
import numpy as np
from math import floor

def calc5thOrderInterp(t,tsol,ysol,fsol,gsol):

    """This function interpolates the results the function value using a quintic
     Hermite spline at t using a vector of function values (ysol), the first 
     derivative (fsol) and second derivative values (gsol) evaluated at (tsol)
    
     @param t   : the argument you want to evaluate y(t) at
     @param tsol: monotonically increasing function argument values
     @param ysol: function values at tsol
     @param fsol: first derivative of function values at tsol
     @param gsol: second derivative of function values at tsol
    
     @returns [y f]: y and f interpolated at t"""
    
    rowMax=tsol.shape[0] 
    
    assert t >= tsol[0] and t <= tsol[rowMax-1],'Error: t not in domain of tsol'

    assert tsol.shape[1] == ysol.shape[1] and tsol.shape[1] == fsol.shape[1],'Error: tsol, ysol, and fsol not the same lengths - they should be'
       
    #Get to the correct subinterval of t
    idx = floor(rowMax/2) #where tsol(idx) <= t
    delta = floor(rowMax/4)
    
    while delta > 4:
       #if t < tsol(idx):
       if t < tsol[idx-1][0]:
           idx = idx - delta
       else:
           idx = idx + delta

       delta = floor(delta/2)

    
    #if t < tsol(idx):
    if t < tsol[idx-1][0]:
       delta = -1
    else:
       delta = 1

    
    #while t >= tsol(idx) and t <= tsol(idx+1) :
    while not (t >= tsol[idx-1][0] and t <= tsol[idx][0]) :
        idx = idx + delta
    
    #Now we have the correct subinterval
    
    idx0 = idx
    idx1 = idx0+1
    
    dtdu = (tsol[idx1-1][0]-tsol[idx0-1][0])
    dudt = 1/dtdu
    
    y0 = ysol[idx0-1][0]
    f0 = fsol[idx0-1][0]*dtdu
    g0 = gsol[idx0-1][0]*dtdu*dtdu
    y1 = ysol[idx1-1][0]
    f1 = fsol[idx1-1][0]*dtdu
    g1 = gsol[idx1-1][0]*dtdu*dtdu
    
    a0 = y0
    a1 = f0
    a2 = 0.5*g0
    
    a345Minv = np.array([[10.0000,-4.0000,0.5000],[-15.0000,7.0000,-1.0000],[6.0000,-3.0000,0.5000]])
    a345RHS = np.array([y1 - (a0+a1+a2),f1 - (a1 + 2*a2),g1 - (2*a2)])
               
    a345 = a345Minv*a345RHS
    a345 =np.dot(a345Minv,a345RHS)

    a3 = a345[0]
    a4 = a345[1]
    a5 = a345[2]
    
    u1 = (t-tsol[idx0-1]) * dudt
    u2 = u1*u1
    u3 = u2*u1
    u4 = u3*u1
    u5 = u4*u1
    
    y =  a0  +  a1*u1 +  a2*u2  +  a3*u3  +  a4*u4  +  a5*u5
    
    return y
