# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calcBezierYFcnXCurveSampleVector.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcBezierYFcnXDerivative import calcBezierYFcnXDerivative
from ComplementaryFonc import minmin,maxmax
import numpy as np

def calcBezierYFcnXCurveSampleVector(curveParams, npts):
    
    """This function evaluates a Bezier spline curve (x(u), y(u)) across its
     entire domain, and and slightly beyond to get values in the extrapolated
     region. The spline is evaluated at its value, and first 3 derivatives. In
     addition if the curve has its integral function defined, then the integral
     is also evaluated.
    
     @param curveParams : a structure with the fields
               .xpts=curveParams[0]
               .ypts=curveParams[1]
               .integral=curveParams[6]
    
       xpts: n x m matrix of Bezier control points defining x(t). 
    
                  n x m matrix 
                  n: Bezier curve order + 1 
                  m: # spline sections
                    
                  Each column defines the control points used for a Bezier 
                  section. Thus a quintic Bezier spline with 3 sections will 
                  have a 6 x 3 matrix for xpts
                  
       ypts: n x m matrix of Bezier control points defining y(t)
    
       integral: if the integral curve has not been numerically computed, 
       then this field will be empty. Otherwise it will have fields of

                      .xptsN : points the integral is evaluated at  
                      .yptsN : numerical value of the integral
                      .y1ptsN: '' integral's first derivative.
                      .y2ptsN: '' integral's 2nd derivative.

     @param npts: the number of samples to use across the curve domain             
     @return curveValues, a struct with the fields"""

    xmin  = minmin(curveParams[0])
    xmax  = maxmax(curveParams[0])
    delta = xmax-xmin
    xmin  = xmin-delta/5
    xmax  = xmax+delta/5
    
    ymin  = minmin(curveParams[1])
    ymax  = maxmax(curveParams[1])
    delta = ymax-ymin
    ymin  = ymin-delta/5
    ymax  = ymax+delta/5
    
    x=np.array([(xmin+(i*(xmax-xmin))/npts) for i in range(npts+1)])
    
    y=np.zeros(shape=x.shape)
    dydx = np.zeros(shape=x.shape)
    d2ydx2 = np.zeros(shape=x.shape)
    d3ydx3 = np.zeros(shape=x.shape)
    intYdx = []
    if curveParams[6] != []:
        intYdx = np.zeros(shape=x.shape)

    for k in range(x.shape[0]):
        if curveParams[6] != []:
            intYdx[k] = calcBezierYFcnXDerivative(x[k], curveParams, -1)
 
        y[k]=calcBezierYFcnXDerivative(x[k], curveParams, 0) 
        dydx[k]   = calcBezierYFcnXDerivative(x[k], curveParams, 1)
        d2ydx2[k] = calcBezierYFcnXDerivative(x[k], curveParams, 2)
        d3ydx3[k] = calcBezierYFcnXDerivative(x[k], curveParams, 3)


    #create curveValues
    curveValues=(x,y,dydx,d2ydx2,d3ydx3,intYdx)
    
    return curveValues
