# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createFiberForceVelocityCurve2018.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcQuinticBezierCornerControlPoints import calcQuinticBezierCornerControlPoints
from math import sqrt
import numpy as np
from numpy import transpose

def createFiberForceVelocityCurve2018(fmaxE,dydxE,dydxC,flag_smoothenNonZeroDyDxC,dydxNearE,fvAtHalfVmax,eccCurviness):
    """ This function will generate a C2 continous (continuous to the second 
     derivative) force velocity curve of a single muscle fiber. The main 
     function of this element is to model the amount the force enhancement or 
     attenuation that is associated with contracting at a particular velocity.
    
     @param fmaxE  The normalized maximum force the fiber can generate when 
            is being stretched. This value is reported to range 
            between 1.1 and 1.8 in the literature, though all values
            are above 1.
            
     @param dydxE    The analogous term of dydxC parameter but for the 
            eccentric portion of the force-velocity curve. As with
            the dydxC term, the physiologically accurate value for
            this parameter is 0, though a value of 0 is rarely used
            in muscle models.  If you are using an equilbrium type 
            model this term must be positive and greater than zero 
            so that the fv curve can be inverted. 
            <br /><br />
            Minimum Value: 0
            Maximum Value: dydxC < (fmaxE-1).
            <br /><br />
            As with the dydxC term, 
            the size of this term also affects the stiffness of the 
            integration problem for equilibrium-type muscle models: 
            the closer to zero this term is, the stiffer the model 
            will be (but only when (dlce(t)/dt)/vmax approaches 1.
     
     @param dydxC  The slope of the fv(dlce(t)/dt) curve at the maximum 
            normalized concentric contraction velocity. Although 
            physiologically the value of dydxC at the maximum 
            concentric contracton velocity is by definition 0, a value
            of 0 is often used. If you are using an equilbrium type 
            model this term must be positive and greater than zero so
            that the fv curve can be inverted.
            <br /><br />
            Minimum Value: 0
            Maximum Value: dydxC < 1 
            <br /><br />
    
     @param flag_smoothenNonZeroDyDxC Setting this flag to 1 will linearly
            extrapolate the force velocity curve to match the terminal
            derivative of the Hill force-velocity curve for shortening 
            contraction velocities faster than vmax. 
    
     @param dydxNearE The slope of the force velocity curve as it approaches
             the maximum eccentric (lengthening) contraction velocity.
             <br /><br />
              Minimum Value: > dydxE
              Maximum Value: dydxNearE < (fmaxE-1)
              <br /><br />
     
     @param fvAtHalfVmax: the value of the force velocity curve at half of
            the maximum shortening velocity. Note that this value must be 
            within the limits listed below.
             <br /><br />
              Minimum Value: > 0.05
              Maximum Value: <= 0.45
              <br /><br />
    
     @param eccCurviness     The dimensionless 'curviness' parameter that 
                can vary between 0 (a line) to 1 (a smooth, but 
                sharply bent elbow). This parameter affects only 
                the eccentric side of the fv curve.
     
     
     @param curveName The name of the muscle this curve applies to. This 
              curve name should have the name of the muscle and the
              curve in it (e.g. 'bicep_fiberForceVelocityCurve') 
              sothat if this curve ever causes an exception, a 
              userfriendly error message can be displayed to the
              end user to help them debug their model.
     
     @throws exception unless these conditions are met
     
      -0 <= dydxC < 1
      -dydxC < dydxNearC < 1
      -1 < dydxIso
      -dydxE < (fmaxE-1) 
      -dydxE < dydxNearC < (fmaxE-1)
      -0<= concCurviness <=0
      -0 <= eccCurviness <= 0
     
     @return fiberForceVelocityCurve 
            A structure that the function calcNormalizedMuscleCurveDerivative 
            can use to evaluate the active force length curve value
            or up to the 3rd derivative."""

    vMax = 1
    eps=2.2204e-16                                  
    vMaxC= -vMax #since internally a concentric velocity is negative because the fiber is getting shorter.
    vMaxE=  vMax
    
    #Ensure that the inputs are within a valid range
    assert fvAtHalfVmax <= 0.45,"must be < 0.45"                   
    
    assert fmaxE > 1.0, "fmaxE must be greater than 1"
    
    assert dydxC >= 0.0 and dydxC < 1,"dydxC must be greater than or equal to 0 and less than 1"
    
    assert dydxE >= 0.0 and dydxE < (fmaxE-1),"dydxE must be greater than or equal to 0 and less than fmaxE-1 (%f)"
    
    assert dydxNearE >= dydxE and dydxNearE < (fmaxE-1),"dydxNearE must be greater than or equal to dydxE and less than fmaxE-1 (%f)"
    
    assert eccCurviness <= 1.0 and eccCurviness >= 0, "eccCurviness must be between 0 and 1"
    
    fiso = 1
    w = 0.5*vMaxC
    a = -fvAtHalfVmax*w*fiso / (vMaxC*fvAtHalfVmax - fiso*vMaxC + fiso*w)
    b =  a*vMaxC/fiso
    
    yCheck  = (b*fiso-a*w)/(b+w);
    assert abs(yCheck-fvAtHalfVmax) < sqrt(eps)
    
    w = 0*vMaxC
    dydxIso =(-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w))
    
    w = 0.9*vMaxC
    dydxNearC = (-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w))
    
    assert abs(dydxNearC) > abs(dydxC) and abs(dydxNearC) < abs(1/vMaxC), "dydxNearC must be greater than or equal to 0 and less than 1"

    #Solve for the concCurviness that results in the Bezier curve that minimizes
    #the error between the Bezier curve and Hill's concentric contraction equations
    #at 10 points between omega = 0, and omega = 0.9*vMax 
    
    xNearC =  0.9*vMaxC
    yNearC = (b*fiso-a*w)/(b+w)
    xIso = 0
    yIso = 1.0
    cC = 0.5
    
    scaleCurviness = lambda curviness:0.1 + 0.8*curviness
    cE = scaleCurviness(eccCurviness)
    
    #Compute the concentric control point locations
    xC= vMaxC
    yC = 0
    if(flag_smoothenNonZeroDyDxC==1):
      dydxC = 0.5*dydxNearC
    
    xNearC = 0.9*vMaxC
    yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC)

    xIso = 0
    yIso = 1
    
    concPts1 =calcQuinticBezierCornerControlPoints(xC,yC,dydxC,0,xNearC, yNearC,dydxNearC,0, cC)
    concPts2 =calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 0,xIso,yIso,dydxIso, 0,cC)

    yIsoH = fmaxE - (0.1*vMaxE*dydxE) - (dydxNearE*0.9*vMaxE)
    dydxIsoH = max(5.0*dydxNearE, 5.0*dydxIso)
    xIsoH = (yIsoH-1.0)/dydxIsoH
    
    xIsoHm = 2*0.5*xIsoH ; yIsoHm  = 1 + (yIsoH-1)*0.5
    
    eccPts1 =calcQuinticBezierCornerControlPoints(xIso,yIso,  dydxIso, 0,xIsoHm, yIsoHm, dydxIsoH, 0, cE)
    
    xE= vMaxE ; yE= fmaxE
    
    xNearE = 0.9*vMaxE ; dydxNearE = dydxNearE/vMaxE
    yNearE = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE)
    
    
    eccPts2 = calcQuinticBezierCornerControlPoints(xIsoHm,yIsoHm,dydxIsoH, 0,xNearE, yNearE, dydxNearE, 0,cE)

    eccPts3 =calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 0,xE,yE,dydxE,0,cE)
    
    xpts =np.array([concPts1[:,0] ,concPts2[:,0] ,eccPts1[:,0], eccPts2[:,0], eccPts3[:,0]])
    ypts = np.array([concPts1[:,1], concPts2[:,1] ,eccPts1[:,1] ,eccPts2[:,1], eccPts3[:,1]])
    
    
    #create fiberForceVelocityCurve
    xpts=transpose(xpts)
    xEnd = [xC, xE]
    yEnd = [yC, yE]
    dydxEnd  = [dydxC, dydxE]
    d2ydx2End= [0,0]
    integral = []
    
    fiberForceVelocityCurve=(xpts,ypts,xEnd,yEnd,dydxEnd,d2ydx2End,integral)
    return fiberForceVelocityCurve
