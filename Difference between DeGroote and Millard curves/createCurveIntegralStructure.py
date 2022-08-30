# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createCurveIntegralStructure.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from calcBezierYFcnXDerivative import calcBezierYFcnXDerivative
from ComplementaryFonc import minmin,maxmax
from scipy.integrate import solve_ivp
import numpy as np

def createCurveIntegralStructure(curveParams, npts, tol, xScaling):

    """Numerically evaluates the integral at npts over the domain of the 2D 
     Bezier curve defined by the matrix of x control points and y control.
     The first and second derivative of the curve integral are also evaluated
     at these points so that it is possible to interpolate the curve using a
     quintic hermine spline.
    
     @param x: value to 
     @param xpts: n x m matrix of Bezier control points defining x(t). 
    
                  n x m matrix 
                  n: Bezier curve order + 1 
                  m: # spline sections
                    
                  Each column defines the control points used for a Bezier 
                  section. Thus a quintic Bezier spline with 3 sections will 
                  have a 6 x 3 matrix for xpts
                  
     @param ypts: n x m matrix of Bezier control points defining y(t)
    
     @param npts: number of intermediate points between xmin and xmax to
                  evaluate the integeral.
     @param tol: relative and absolute tolerance on the integral.
    
     @return integralStruct: A structure containing the numerically calculated
                             integral, its first and second derivative so that 
                             the integral curve can be interpolated using a 
                             quintic Hermite spline.
    
               Has fields of
                          .xptsN : points the integral is evaluated at  
                          .yptsN : numerical value of the integral
                          .y1ptsN: '' integral's first derivative.
                          .y2ptsN: '' integral's 2nd derivative."""

    fcn = lambda t,x: calcBezierYFcnXDerivative(t, curveParams, 0)
    x=curveParams[0]
    
    xmin = minmin(x)
    xmax = maxmax(x)
    
    xv=np.array([xmin+i*(xmax-xmin)/(npts-1) for i in range(npts)] ) 
    
    sol=solve_ivp(fcn,[xmin,xmax],[0],rtol=tol,atol=tol,dense_output=True)
    ye=sol.sol(xv) 
    
    xptsN=np.resize(xv,new_shape=(len(xv),1))
    yptsN=np.resize(ye,new_shape=(np.shape(ye)[1],1))
    
    y1ptsN = np.zeros(shape=yptsN.shape)
    y2ptsN = np.zeros(shape=yptsN.shape)

    for i in range(yptsN.shape[1]):
       y1ptsN[0][i] =  calcBezierYFcnXDerivative(xptsN[i], curveParams, 0) 
       y2ptsN[0][i] =  calcBezierYFcnXDerivative(xptsN[i], curveParams, 1) 
    
    integralStruct=(xptsN,yptsN,y1ptsN,y2ptsN,xScaling)
   
    return integralStruct

