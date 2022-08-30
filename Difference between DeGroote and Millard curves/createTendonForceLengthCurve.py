# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createTendonForceLengthCurve.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcQuinticBezierCornerControlPoints import calcQuinticBezierCornerControlPoints
from createCurveIntegralStructure import createCurveIntegralStructure
import numpy as np
from numpy import transpose

def createTendonForceLengthCurve( eIso, kIso,fToe, curviness,computeIntegral):
    """%Will generate a C2 continous (continuous to the second derivative) 
    curve in a MuscleFunctionObject object that fits a tendon's tensile 
    force length curve. 
    
    
    
    @param eIso   The tendon strain at which the tendon develops 1 unit
    			of normalized force (1 maximum isometric force). Note that 
    			the'1' is left off. Thus an e0 of 0.04 means that the tendon 
    			will develop an 1 normalized force unit when it is strained 
    			by 4% of its resting length, at a normalized length of 
    			1.04
    
    @param kIso    The normalized stiffness (or slope) of the tendon
    				curve when the tendon is strained by e0 
    				(or has a length of 1+e0) under a load of 1 maximum
    				isometric unit of force.        
    
    @param fToe    The normalized force at which the tendon smoothly
    			   transitions from the curved low stiffness region to 
    			   the linear stiffness region.
    
    @param curviness    The dimensionless 'curviness' parameter that 
    					can vary between 0 (a line) to 1 (a smooth, but 
    					sharply bent elbow)
    
    @param computeIntegral  If this is true, the integral for this curve
    						is numerically calculated and splined. If false, 
    						this integral is not computed, and a call to 
    						.calcIntegral will throw an exception
    
     @param curveName The name of the muscle this curve applies to. This 
    				  curve name should have the name of the muscle and the
    				  curve in it (e.g. 'bicep_tendonForceLengthCurve') 
    				  sothat if this curve ever causes an exception, a 
    				  userfriendly error message can be displayed to the
    				  end user to help them debug their model.
    
    @throws SimTK::Exception unless the following conditions are met:
    	-0 < fToe < 1
    	-e0 > 0
    	-kiso > 1/e0
    	-0 <= curviness <= 1
    
    @return tendonForceLengthCurve"""

    #Check the input arguments
    assert eIso>0 ,"eIso must be greater than 0"
    
    assert fToe>0 and fToe < 1 ,"fToe must be greater than 0 and less than 1, but %f was entered"
    
    assert kIso > (1/eIso) , " kIso must be greater than 1/eIso, (%f), but kIso (%f) was entered"
    
    assert curviness>=0 and curviness <= 1 ,"curviness must be between 0.0 and 1.0, but %f was entered"
   
    #Translate the user parameters to quintic Bezier points
    scaleCurviness= lambda curviness:0.1 + 0.8*curviness
    c = scaleCurviness(curviness)
    x0 = 1.0 
    y0 = 0 
    dydx0 = 0
    
    xIso = 1.0 + eIso 
    yIso = 1 ; dydxIso = kIso
    
    #Location where the curved section becomes linear
    yToe = fToe ;
    xToe = (yToe-1)/kIso + xIso
    
    
    #To limit the 2nd derivative of the toe region the line it tends to
    #has to intersect the x axis to the right of the origin
    xFoot = 1.0+(xToe-1.0)/10.0
    yFoot = 0
    
    #Compute the location of the corner formed by the average slope of the
    #toe and the slope of the linear section
    yToeMid = yToe*0.5 ; xToeMid = (yToeMid-yIso)/kIso + xIso
    dydxToeMid = (yToeMid-yFoot)/(xToeMid-xFoot)
    
    #Compute the location of the control point to the left of the corner
    xToeCtrl = xFoot + 0.5*(xToeMid-xFoot); 
    yToeCtrl = yFoot + dydxToeMid*(xToeCtrl-xFoot)

    #Compute the Quintic Bezier control points
    p0 = calcQuinticBezierCornerControlPoints(x0,y0,dydx0, 0,xToeCtrl,yToeCtrl,dydxToeMid, 0, c)
    p1 = calcQuinticBezierCornerControlPoints(xToeCtrl, yToeCtrl, dydxToeMid, 0,xToe,yToe,dydxIso, 0, c)
    xpts = np.array([p0[:,0],p1[:,0]])
    ypts = np.array([p0[:,1],p1[:,1]])
    
    #create tendonForceLengthCurve
    xpts=transpose(xpts)
    xEnd= [x0, xToe]
    yEnd= [y0, yToe]
    dydxEnd= [dydx0, dydxIso]
    d2ydx2End= [0, 0]
    integral = []
    
    tendonForceLengthCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    
    if(computeIntegral == 1):
        xScaling = eIso
        integral =createCurveIntegralStructure(tendonForceLengthCurve,1000,1e-12,xScaling);  
        
    tendonForceLengthCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    
    return tendonForceLengthCurve

