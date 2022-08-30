# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createFiberActiveForceLengthCurve.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcQuinticBezierCornerControlPoints import calcQuinticBezierCornerControlPoints
from math import sqrt
import numpy as np
from numpy import transpose
        
def createFiberActiveForceLengthCurve(lce0, lce1, lce2, lce3,minActiveForceLengthValue, plateauSlope,curviness,computeIntegral):  
    """       This is a function that will produce a C2 (continuous to the second
           derivative) active force length curve.
    
    
           @param lce0   Normalized fiber length at the left-most shoulder of the 
                         active force-length curve. The value of the active force
                         length curve for lce < lce0 will be equal to the value
                         set in shoulderVal. Normally lce0 is approximately 0.5
           
           @param lce1   Normalized fiber length at the transition point between 
                         the ascending limb and the plateau region of the active 
                         force length curve.
           
           @param lce2   Normalized fiber length at the maximum active force length
                         curve value of 1. Normally lce2 is by definition 1.
           
           @param lce3   Normalized fiber length of the at the right most shoulder
                         of the active-force length curve. The value of the active
                         force length curve for lce > lce2 will be equal to the 
                         value of shoulderVal. Normally lce3 is approximately 1.5
    
           @param minActiveForceLengthValue
                                 The minimum value of the active force length 
                                 curve. A physiological non-equibrium muscle model
                                 would have this value set to 0. An equilibrium 
                                 muscle model would have a non-zero lower bound on 
                                 this value of 0.1 typically. shoulderVal must be 
                                 greater than, or equal to 0.
                               
           @param plateauSlope   The slope of the plateau of the active force
                                 length curve between lce1 and lce2. This parameter
                                 can vary depending on the muscle model, but a 
                                 value of 0.8616 is a good place to start.
    
           @param curviness  The dimensionless 'curviness' parameter that 
                             can vary between 0 (a line) to 1 (a smooth, but 
                             sharply bent elbow). A value of 0 will yield an active 
                             force length curve that is composed of slightly curved 
                             line segments. A value of 1 will yield an active force
                             length curve that is smoothly rounded.
    
           @param computeIntegral If this is true, the integral for this curve
                                  is numerically calculated and splined. If false, 
                                  this integral is not computed, and a call to 
                                  .calcIntegral will throw an exception
    
           @param muscleName The name of the muscle this curve applies to. This 
                             curve name should have the name of the muscle and the
                             curve in it (e.g. "bicep_fiberActiveForceLengthCurve") 
                             sothat if this curve ever causes an exception, a 
                             userfriendly error message can be displayed to the
                             end user to help them debug their model.
    
           @throws SimTK::Exception if these conditions aren't met
               -0 < lce0 < lce1 < lce2 < lce3 
               -shoulderVal >= 0
               -0 <= plateauSlope < (1/(lce3-lce2))
               -0 <= curviness <= 1
    
           @return activeForceLengthCurve 
                   A structure that the function calcNormalizedMuscleCurveDerivative can use  
                   to evaluate the active force length curve value or up to the 3rd derivative."""


    x0=lce0
    x1=lce1
    x2=lce2 
    x3=lce3
    ylow=minActiveForceLengthValue
    dydx=plateauSlope
    
    #------------
    #Check inputs
    #------------
    eps=2.2204e-16 
    rootEPS=sqrt(eps)
    assert x0>=0 and x1>x0+rootEPS and x2>x1+rootEPS and x3>x2+rootEPS,"This must be true: 0 < lce0 < lce1 < lce2 < lce3"
    
    assert ylow>= 0, "shoulderVal must be greater than, or equal to 0"
    
    dydxUpperBound=(1-ylow)/(x2-x1)
    assert dydx >= 0 and dydx < dydxUpperBound , "plateauSlope must be greater than 0 and less than"

    assert curviness >= 0 and curviness <= 1 , "curviness must be between 0 and 1"

    #----------------
    #Create the curve
    #----------------
    #Translate the users parameters into Bezier curves 
    scaleCurviness= lambda curviness:0.1 + 0.8*curviness
    c = scaleCurviness(curviness)
    #Calculate the location of the shoulder
    xDelta = 0.05*x2 #half the width of the sarcomere 0.0259, 
                       #but TM.Winter's data has a wider shoulder than this
    xs=(x2-xDelta)

    #Calculate the intermediate points located on the ascending limb
    y0=0
    dydx0=0
    y1=1-dydx*(xs-x1)
    dydx01= 1.25*(y1-y0)/(x1-x0)
    x01=x0 + 0.5*(x1-x0)
    y01= y0 + 0.5*(y1-y0)
    #Calculate the intermediate points of the shallow ascending plateau
    x1s=x1 + 0.5*(xs-x1)
    y1s=y1 + 0.5*(1-y1)
    dydx1s= dydx
    
    y2 = 1
    dydx2 = 0 
    
    #Descending limb
    y3 = 0 
    dydx3 = 0 #x3 entered
    
    x23 = (x2+xDelta) + 0.5*(x3-(x2+xDelta))
    y23 = y2 + 0.5*(y3-y2)
    dydx23  = (y3-y2)/((x3-xDelta)-(x2+xDelta)); 
    
    #Compute the locations of the control points
    p0 = calcQuinticBezierCornerControlPoints(x0,ylow,dydx0,0,x01,y01,dydx01,0,c)
                                         
    p1 = calcQuinticBezierCornerControlPoints(x01,y01,dydx01,0,x1s,y1s,dydx1s,0,c)
                                         
    p2 = calcQuinticBezierCornerControlPoints(x1s,y1s,dydx1s,0,x2,y2,dydx2,0,c)
                                          
    p3 = calcQuinticBezierCornerControlPoints(x2,y2,dydx2,0,x23,y23,dydx23,0,c)
                                         
    p4 = calcQuinticBezierCornerControlPoints(x23,y23,dydx23,0,x3,ylow,dydx3,0,c)
    
    xpts=np.array([p0[:,0], p1[:,0], p2[:,0], p3[:,0] ,p4[:,0]])
    ypts=np.array([p0[:,1], p1[:,1], p2[:,1], p3[:,1] ,p4[:,1]])
    
    #Create the curve structure
    xpts=transpose(xpts)
    xEnd = [x0, x3]
    yEnd = [ylow, ylow]
    dydxEnd = [dydx0, dydx3]
    d2ydx2End = [0, 0]
    integral = []

    activeForceLengthCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    
    return activeForceLengthCurve
