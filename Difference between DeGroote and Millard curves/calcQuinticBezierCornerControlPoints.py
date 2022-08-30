# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calcQuinticBezierCornerControlPoints.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from math import sqrt
import numpy as np

def calcQuinticBezierCornerControlPoints(x0,y0,dydx0,d2ydx20,x1,y1,dydx1,d2ydx21,curviness):
    
    """Note: This is an improved version of the code that is in OpenSim: this
           function allows you to set the 2nd derivative at the end points of
           the Bezier curve section. The code in OpenSim automatically sets
           these endpoints to 0.
    
            Calculates the location of quintic Bezier curve control points to 
            create a C shaped curve.
    
            @param x0       First intercept x location
            @param y0       First intercept y location
            @param dydx0    First intercept slope
            @param x1       Second intercept x location
            @param y1       Second intercept y location
            @param dydx1    Second intercept slope
            @param curviness A parameter that ranges between 0 and 1 to denote a 
                             straight line or a curve
            @throws OpenSim::Exception 
             -If the curviness parameter is less than 0, or greater than 1;
             -If the points and slopes are chosen so that an "S" shaped curve would 
              be produced. This is tested by examining the points (x0,y0) and 
              (x1,y1) together with the intersection (xC,yC) of the lines beginning 
              at these points with slopes of dydx0 and dydx1 form a triangle. If the 
              line segment from (x0,y0) to (x1,y1) is not the longest line segment, 
              an exception is thrown. This is an overly conservative test as it 
              prevents very deep 'V' shapes from being respresented.
    
            @return a SimTK::Matrix of 6 points Matrix(6,2) that correspond to the 
                             X, and Y control points for a quintic Bezier curve that
                             has the above properties"""
    
    
    #Check inputs
    xyPts = np.zeros(shape=(6,2))
    assert curviness >= 0 and curviness <= 1 , 'Error: curviness must be [0 1]'      

    xC = 0 
    yC = 0
    eps=2.2204e-16 
    rootEPS=sqrt(eps)
    
    if(abs(dydx0-dydx1) > rootEPS):
        xC = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);   
    else:
        xC = (x1+x0)/2   
    yC = (xC-x1)*dydx1 + y1

    #Start point 
    xyPts[0,0]=x0 
    xyPts[0,1]=y0
    
    #End point 
    xyPts[5,0]=x1
    xyPts[5,1]=y1
    
    #Original code - leads to 2 localized corners
    xyPts[1,0]=x0 + curviness*(xC-xyPts[0,0])
    xyPts[1,1]=y0 + curviness*(yC-xyPts[0,1])

    dxdu0= 5*(xyPts[1,0]-xyPts[0,0])

    xyPts[2,0] = xyPts[1,0] + 0.5*(xC-xyPts[1,0])
    
    d2xdu20 = 20*(xyPts[2,0] - 2*xyPts[1,0] + xyPts[0,0])
    d2ydu20 = (dxdu0*dxdu0*(d2ydx20) + d2xdu20*(dydx0))
    
    xyPts[2,1] = d2ydu20*(1/20) + 2*xyPts[1,1] - xyPts[0,1] 
    
    xyPts[4,0] = xyPts[5,0] + curviness*(xC-xyPts[5,0])
    
    xyPts[4,1] = xyPts[5,1] + curviness*(yC-xyPts[5,1])
    
    dxdu1 = 5*(xyPts[5,0]-xyPts[4,0])
    
    xyPts[3,0] = xyPts[4,0] + 0.5*(xC-xyPts[4,0])
    
    d2xdu21 = 20*(xyPts[3,0] - 2*xyPts[4,0] + xyPts[5,0])

    d2ydu21 = (dxdu1*dxdu1*(d2ydx21) + d2xdu21*(dydx1))

    xyPts[3,1] = d2ydu21*(1/20) + 2*xyPts[4,1] - xyPts[5,1] 
    
    return xyPts