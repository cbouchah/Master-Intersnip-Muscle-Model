# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createInverseBezierCurve.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from numpy import transpose

def createInverseBezierCurve(curve):
    
    Inv_xpts = transpose(curve[1])
    Inv_ypts = curve[0]
    Inv_xEnd = curve[3]
    Inv_yEnd = curve[2]
    
    col  = curve[0].shape[1] 
    row  = curve[0].shape[0]
    n    = row-1 #polynomial order

    dx1 = n*(Inv_xpts[1:row,:]-Inv_xpts[0:(row-1),:])
    dy1 = n*(Inv_ypts[1:row,:]-Inv_ypts[0:(row-1),:])

    dx2 = (n-1)*(dx1[1:(row-1),:]-dx1[0:(row-2),:])
    dy2 = (n-1)*(dy1[1:(row-1),:]-dy1[0:(row-2),:])
    
    dydx0 = dy1[0,0]/dx1[0,0]
    dydx1 = dy1[row-2,col-1]/dx1[row-2,col-1]
    
    d2ydx20 = 0
    d2ydx21 = 0

    if(curve[5][0] != 0):
        d2ydx20 = (dy2[0,0]*dx1[0,0]- dy1[0,0]*dx2[0,0])/( dx1[0,0]**2 )

    if(curve[5][1]!= 0):
        d2ydx21 = (dy2[row-3,col-1]*dx1[row-2,col-1]- dy1[row-2,col-1]*dx2[row-3,col-1])/( dx1[row-2,col-1]**2 )

    #create curveInv
    Inv_ypts = transpose(curve[0])
    Inv_dydxEnd   = [  dydx0,  dydx1]
    Inv_d2ydx2End = [d2ydx20, d2ydx21]

    Inv_integral = []
    
    curveInv=(Inv_xpts,Inv_ypts,Inv_xEnd,Inv_yEnd,Inv_dydxEnd,Inv_d2ydx2End,Inv_integral)
    
    return curveInv
