# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createFiberForceLengthCurve.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcQuinticBezierCornerControlPoints import calcQuinticBezierCornerControlPoints
from createCurveIntegralStructure import createCurveIntegralStructure
import numpy as np
from numpy import transpose

def createFiberForceLengthCurve(eZero, eIso,kLow, kIso, curviness,computeIntegral):
    """This function will generate a C2 continuous curve that fits a fiber's 
    tensile force length curve.
    
    @param eZero The fiber strain at which the fiber begins to develop force.
    			 Thus an e0 of 0.0 means that the fiber will start to develop
    			 passive force when it has a normalized length of 1.0. Note
    			 that e0 can be postive or negative.
    
    @param eIso The fiber strain at which the fiber develops 1 unit of 
    			normalized force (1 maximum isometric force). Note that the 
    			'1' is left off. Thus an e0 of 0.6 means that the fiber 
    			will develop an 1 normalized force unit when it is strained 
    			by 60% of its resting length, or to a normalized length of 
    			1.6
    
    @param kLow   The normalized stiffness (or slope) of the fiber curve 
    			  close to the location where the force-length curve 
    			  approaches a normalized force of 0. This is usually 
    			  chosen to be a small, but non-zero fraction of kIso 
    			  (kLow = 0.025 kIso is typical).
    
    @param kIso   The normalized stiffness (or slope) of the fiber curve 
    			  when the fiber is strained by eIso (or has a length of 
    			  1+eIso) under a load of 1 maximum isometric unit of force.
    
    
    @param curviness    The dimensionless 'curviness' parameter that 
    					can vary between 0 (a line) to 1 (a smooth, but 
    					sharply bent elbow)
    
    @param computeIntegral  If this is true, the integral for this curve
    						is numerically calculated and splined. If false, 
    						this integral is not computed, and a call to 
    						.calcIntegral will throw an exception
    
     @param curveName The name of the muscle this curve applies to. This 
    				  curve name should have the name of the muscle and the
    				  curve in it (e.g. "bicep_fiberForceLengthCurve") 
    				  so that if this curve ever causes an exception, a 
    				  userfriendly error message can be displayed to the
    				  end user to help them debug their model.
    
    @throws exception unless the following conditions are met
    	-eIso > eZero            
    	-kIso > 1/(eIso-eZero)
    	-0 < kLow < kIso
    	-0 <= curviness <= 1
    
    @return fiberForceLengthCurve 
           A structure that the function calcNormalizedMuscleCurveDerivative 
           can use to evaluate the active force length curve value
           or up to the 3rd derivative."""

    #Check the input arguments
    assert eIso > eZero ,"The following must hold: eIso  > eZero"
    
    assert kIso > 1.0/(eIso-eZero) , "kIso must be greater than 1/(eIso-eZero)"
    
    assert kLow > 0.0 and kLow < 1/(eIso-eZero) , "kLow must be greater than 0 and less than or equal to 1"
    
    assert curviness>=0 and curviness <= 1, "curviness must be between 0.0 and 1.0"
    
    #Translate the user parameters to quintic Bezier points
    scaleCurviness= lambda curviness:0.1 + 0.8*curviness
    c = scaleCurviness(curviness)
    xZero = 1+eZero
    yZero = 0
    xIso = 1 + eIso
    yIso = 1
    
    deltaX = min(0.1*(1.0/kIso), 0.1*(xIso-xZero))
    
    xLow= xZero + deltaX 
    xfoot= xZero + 0.5*(xLow-xZero)
    yfoot = 0
    yLow = yfoot + kLow*(xLow-xfoot)
    
    #Compute the Quintic Bezier control points
    p0 = calcQuinticBezierCornerControlPoints(xZero, yZero,0, 0,xLow, yLow, kLow, 0,c)
    
    p1 = calcQuinticBezierCornerControlPoints(xLow, yLow, kLow, 0,xIso, yIso, kIso, 0, c)
    xpts = np.array([p0[:,0], p1[:,0]])
    ypts = np.array([p0[:,1], p1[:,1]])
    
    #Create the curve structure

    xpts=transpose(xpts)
    xEnd = [xZero, xIso]
    yEnd = [yZero, yIso]
    dydxEnd= [0, kIso]
    d2ydx2End = [0, 0]
    integral = []

    fiberForceLengthCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    
    if(computeIntegral == 1):
        xScaling = eIso
        integral = createCurveIntegralStructure(fiberForceLengthCurve,1000,1e-12,xScaling)   
       
    fiberForceLengthCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    
    return fiberForceLengthCurve
