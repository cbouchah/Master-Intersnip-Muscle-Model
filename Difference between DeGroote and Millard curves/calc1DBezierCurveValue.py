# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calc1DBezierCurveValue.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
import numpy as np

def  calc1DBezierCurveValue(u, pV):

    """ This function implements De Casteljau's recursive algorithm for
     evaluating an nth order 1D Bezier curve of the form
    
      B(u) = sum_{i=0}^n  [(n choose i)(1-u)^{n-1} u^i] p_i
    
     where 
      u: argument of the curve
      n: order of the curve
      p_i: value of the ith control point
    
     For an n th order curve with (n+1) points this algorithm requires n! 
     subtractions, multiplications, and additions to terminate, and a stack 
     that is n! deep. Although this algorithm is very general faster results 
     can be obtained using optimized code if the order of the Bezier curve is 
     known ahead of time using optimized code.
    
     @params u : [0,1] the argument of the Bezier curve
     @params pV: vector of control points
    
     @returns f: the value of the Bezier curve evaluated at u"""

    bV0=np.zeros(len(pV)-1) 
    
    for i in range(len(bV0)):
        bV0[i] = (pV[i+1]-pV[i])*u + pV[i]

    if len(bV0)==1:
        f  = bV0[0]
    else:
        f = calc1DBezierCurveValue(u,bV0)
        
    return f