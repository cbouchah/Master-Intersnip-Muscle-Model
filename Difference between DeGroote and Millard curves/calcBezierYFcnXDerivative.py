# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calcBezierYFcnXDerivative.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
import numpy as np
from numpy import random
from calcIndex import calcIndex
from calc1DBezierCurveValue import calc1DBezierCurveValue
from calc5thOrderInterp import calc5thOrderInterp
from ComplementaryFonc import minmin, maxmax

def calcBezierYFcnXDerivative(x, curveParams, der):

    """ This function takes the control points for the Bezier curves
     x(t) and f(t) and treats it as a function:
    
            val = y(x)
    
     To do this we assume that
     1. x(t) is a monotonic increasing function
     2. d/dt x(t) > 0
    
     In addition, to maintain continuity we assume that the second derivative
     of the curve goes to 0 at the end points so that we can linearly
     extrapolate the curve outside of [xmin, xmax]
    
     @param x: value to 
    
     @param curveParams=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2E,integral)
    
       xpts: n x m matrix of Bezier control points defining x(t). 
    
                  n x m matrix 
                  n: Bezier curve order + 1 
                  m: # spline sections
                    
                  Each column defines the control points used for a Bezier 
                  section. Thus a quintic Bezier spline with 3 sections will 
                  have a 6 x 3 matrix for xpts
                  
       ypts: n x m matrix of Bezier control points defining y(t)
    
       integral: if the integral curve has not been numerically computed, 
                 then this field will be empty. Otherwise it will have fields
                of
    
                          .xptsN : points the integral is evaluated at  
                          .yptsN : numerical value of the integral
                          .y1ptsN: '' integral's first derivative.
                          .y2ptsN: '' integral's 2nd derivative.
    
     @param der: integer defining the order of the desired derivative. 
                 The value of der is restricted to:
                 
                 0: y = f(x)
                 1: dy/dx  
                 2: d^2y/dx^2 
                 3: d^3y/dx^3
    
    @returns y : d^n/dx^n y(x), the n^th derivative of the function y(x)"""

    val=None
    
    assert der >= -1 and der <= 3,'der must be within [0,3]'
    
    #nrow = size(xpts,1)
    #col = size(xpts,2)
    nrow = curveParams[0].shape[0]
    #ncol = curveParams[0].shape[1]
    xmin = minmin(curveParams[0])
    xmax = maxmax(curveParams[0])
    
    if x < xmin or x > xmax:
        if x <= xmin:
           idxEnd = 0 
        else:
           idxEnd = 1 

        #evaluate the desired derivative of y at x

        if der==-1:
            assert curveParams[6]!=[],'Integral function has not been computed for this curve'
            
            if x <= xmin:
                val = 0
            else:

                y0  = curveParams[6][1][len(np.array(curveParams[6][0]))-1][0]
                f1  = curveParams[6][1][len(np.array(curveParams[6][0]))-1][0]
                f2  = curveParams[6][1][len(np.array(curveParams[6][0]))-1][0]
                x0  = curveParams[6][1][len(np.array(curveParams[6][0]))-1][0]
                   
                val = y0 + f1*(x-x0) + (0.5*f2)*(x-x0)**2
       
        elif der==0:
            x0   = curveParams[2][idxEnd] 
            y0   = curveParams[3][idxEnd]
            dydx = curveParams[4][idxEnd]
                            
            val = dydx*(x-x0) + y0
                 
        elif der==1:
            dydx = curveParams[4][idxEnd]
            val  = dydx
            
        else:
            val = 0 
      
    else:
        #Find the spline section that x is in 
        moduloRange = [] 
        eps=2.2204e-16
        tol= eps        
        #options=[moduloRange,tol]
        xInterval= np.array([curveParams[0][0,:],curveParams[0][nrow-1,:]])
        col=calcIndex(x,xInterval,moduloRange,tol)
                                   
        #Extract the vector of control points, and calculate
        #the control points that define the derivative curve

        xV  = curveParams[0][:,col]       
        x1V = np.diff(xV) *(nrow-1)
        
        #Find the value of u that corresponds to the desired value of x

        u = (x-curveParams[0][0,col]) / (curveParams[0][nrow-1, col]-curveParams[0][0,col])
        iter= 1
        iterMax = 100
        tol = eps*10
        err = tol*10
        
        while iter < iterMax and abs(err) > tol:
           err = calc1DBezierCurveValue(u, xV) - x
           derr = calc1DBezierCurveValue(u, x1V)
           
           if abs(err) > tol and abs(derr) > eps :  
               du = -err/derr
               u  = u + du
              
               #For very nonlinear curves Newton's method can
               #become pathological. If u is outside [0,1] we
               #kick it back into the interval by some random small amount.
               
               if u < 0 or u > 1.0:
                    if u < 0:
                       u = 0 
                       u = u + random.uniform(0,1,(1,1))*0.1 
                    else:
                       u = 1 
                       u = u - random.uniform(0,1,(1,1))*0.1 
         
           iter = iter+1
            
        #Evaluate the desired derivative of y   

        if der==-1:
            assert curveParams[6]!=[],'Integral function has not been computed for this curve'
               
            y0 = calc5thOrderInterp(x,np.array(curveParams[6][0]),np.array(curveParams[6][1]),np.array(curveParams[6][2]),np.array(curveParams[6][3]))
            val = y0 
                  
        elif der==0:
            yV  = curveParams[1][col,:]
            val = calc1DBezierCurveValue(u, yV)
        elif der==1:
            yV  = curveParams[1][col,:]
            y1V =  np.diff(yV) *(nrow-1)
            x1  = calc1DBezierCurveValue(u, x1V)
            y1  = calc1DBezierCurveValue(u, y1V) 
            val = y1/x1
            
        elif der==2:
            x2V =  np.diff(x1V)*(nrow-2) 
                
            yV  =  curveParams[1][col,:]
            y1V =  np.diff(yV) *(nrow-1) 
            y2V =  np.diff(y1V)*(nrow-2) 
        
            x1  = calc1DBezierCurveValue(u, x1V)
            y1  = calc1DBezierCurveValue(u, y1V)
            x2  = calc1DBezierCurveValue(u, x2V)
            y2  = calc1DBezierCurveValue(u, y2V)
                
            t1 = 1/x1
            t3 = x1*x1
                
            val = (y2 * t1 - y1 / t3 * x2) * t1
            
        else: 

            x2V = np.diff(x1V)*(nrow-2) 
            x3V = np.diff(x2V)*(nrow-3) 
                
            yV  =  curveParams[1][col,:] 
            y1V =  np.diff(yV) *(nrow-1) 
            y2V =  np.diff(y1V)*(nrow-2) 
            y3V =  np.diff(y2V)*(nrow-3) 
                
            x1  = calc1DBezierCurveValue(u, x1V)
            y1  = calc1DBezierCurveValue(u, y1V)
                
            x2  = calc1DBezierCurveValue(u, x2V)
            y2  = calc1DBezierCurveValue(u, y2V)
                
            x3  = calc1DBezierCurveValue(u, x3V)
            y3  = calc1DBezierCurveValue(u, y3V)
            
            t1 = 1 / x1
            t3 = x1*x1
            t4 = 1 / t3
            t11 = x2*x2
            t14 = y1 * t4
                
            val = ((y3*t1 - 2*y2*t4*x2+ 2*y1/t3/x1 * t11 - t14 * x3) * t1-(y2*t1 - t14*x2)*t4*x2) * t1      

    return val

def calcBezierYFcnXDerivative2(x, curveParams):

    """ This function takes the control points for the Bezier curves
     x(t) and f(t) and treats it as a function:
    
            val = y(x)
    
     To do this we assume that
     1. x(t) is a monotonic increasing function
     2. d/dt x(t) > 0
    
     In addition, to maintain continuity we assume that the second derivative
     of the curve goes to 0 at the end points so that we can linearly
     extrapolate the curve outside of [xmin, xmax]
    
     @param x: value to 
    
     @param curveParams=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2E,integral)
    
       xpts: n x m matrix of Bezier control points defining x(t). 
    
                  n x m matrix 
                  n: Bezier curve order + 1 
                  m: # spline sections
                    
                  Each column defines the control points used for a Bezier 
                  section. Thus a quintic Bezier spline with 3 sections will 
                  have a 6 x 3 matrix for xpts
                  
       ypts: n x m matrix of Bezier control points defining y(t)
    
       integral: if the integral curve has not been numerically computed, 
                 then this field will be empty. Otherwise it will have fields
                of
    
                          .xptsN : points the integral is evaluated at  
                          .yptsN : numerical value of the integral
                          .y1ptsN: '' integral's first derivative.
                          .y2ptsN: '' integral's 2nd derivative.
    
     @param der: integer defining the order of the desired derivative. 
                 The value of der is restricted to:
                 
                 0: y = f(x)
                 1: dy/dx  
                 2: d^2y/dx^2 
                 3: d^3y/dx^3
    
    @returns y : d^n/dx^n y(x), the n^th derivative of the function y(x)"""

    val=None
    

    nrow = curveParams[0].shape[0]
    xmin = minmin(curveParams[0])
    xmax = maxmax(curveParams[0])
    
    if x < xmin or x > xmax:
        if x <= xmin:
           idxEnd = 0 
        else:
           idxEnd = 1 

        #evaluate the desired derivative of y at x

        x0   = curveParams[2][idxEnd] 
        y0   = curveParams[3][idxEnd]
        dydx = curveParams[4][idxEnd]
                        
        val = dydx*(x-x0) + y0
       
      
    else:
        #Find the spline section that x is in 
        moduloRange = [] 
        eps=2.2204e-16
        tol= eps        
        #options=[moduloRange,tol]
        xInterval= np.array([curveParams[0][0,:],curveParams[0][nrow-1,:]])
        col=calcIndex(x,xInterval,moduloRange,tol)
                                   
        #Extract the vector of control points, and calculate
        #the control points that define the derivative curve

        xV  = curveParams[0][:,col]       
        x1V = np.diff(xV) *(nrow-1)
        
        #Find the value of u that corresponds to the desired value of x

        u = (x-curveParams[0][0,col]) / (curveParams[0][nrow-1, col]-curveParams[0][0,col])
        iter= 1
        iterMax = 100
        tol = eps*10
        err = tol*10
        
        while iter < iterMax and abs(err) > tol:
           err = calc1DBezierCurveValue(u, xV) - x
           derr = calc1DBezierCurveValue(u, x1V)
           
           if abs(err) > tol and abs(derr) > eps :  
               du = -err/derr
               u  = u + du
              
               #For very nonlinear curves Newton's method can
               #become pathological. If u is outside [0,1] we
               #kick it back into the interval by some random small amount.
               
               if u < 0 or u > 1.0:
                    if u < 0:
                       u = 0 
                       u = u + random.uniform(0,1,(1,1))*0.1 
                    else:
                       u = 1 
                       u = u - random.uniform(0,1,(1,1))*0.1 
         
           iter = iter+1
            
        #Evaluate the desired derivative of y   
        yV  = curveParams[1][col,:]
        val = calc1DBezierCurveValue(u, yV)

    return val

