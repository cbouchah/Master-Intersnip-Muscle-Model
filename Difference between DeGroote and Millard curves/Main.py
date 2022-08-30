# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------------------------------------------------
#DE GROOTE MODEL
#From "Evaluation of Direct Collocation Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem"
#Written by: Friedl De Groote, Allison L. Kinney, Anil V. Rao, Benjamin J. Fregly
#--> All numerical values ​​are available on the Online Supplement
#----------------------------------------------------------------------------------------------------------------------

from Mathematical_expressions_for_muscle_tendon_characteristics_De_Groote import f_a,f_p,f_v,f_t

#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> All numerical values ​​are available on the online supplement:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createDefaultNormalizedMuscleCurves.m 
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from createFiberActiveForceLengthCurve import createFiberActiveForceLengthCurve
from createFiberPassiveForceLengthCurve import createFiberForceLengthCurve
from createFiberForceVelocityCurve2018 import createFiberForceVelocityCurve2018
from createTendonForceLengthCurve import createTendonForceLengthCurve
from createInverseBezierCurve import createInverseBezierCurve

from math import log
import numpy as np

from createDefaultNormalizedMuscleCurves import ModelCurvesComparison

""" Construction of all the curves: model of De Groote, of Millard and comparison
Please choose save=True if you want to save the curves """

npts=100
save=False

#Active force length curve

lce0=0.47-0.0259 
lce1=0.73 
lce2=1.0 
lce3=1.8123
minActiveForceLengthValue=0 
curviness=1.0
computeIntegral=0 
plateauSlope=0.8616

l0_Ma=0.2
lf_Ma=2

l_MaG=np.linspace(l0_Ma,lf_Ma,npts)
FA=[f_a(x) for x in l_MaG]
activeForceLengthCurve_DG=[l_MaG,f_a(l_MaG),l0_Ma,lf_Ma,min(FA),max(FA)]
activeForceLengthCurve_M = createFiberActiveForceLengthCurve(lce0,lce1,lce2,lce3,minActiveForceLengthValue,plateauSlope,curviness,computeIntegral)
activeForceLengthParams=["Active Force Length Curve","$f_{act}(l_{MN})$","$l_{MN}$","$f_{act}$","darkred","r","salmon"]

#--------------------------------------
#Hacked Fiber Active Force Length Curve

lce0 = 0.47-0.0259
lce1 = 0.73
lce2 = 1.0
lce3 = 1.8123
minActiveForceLengthValueHack = 0.1 #Here's the hack: minimum value > 0
curviness = 1.0                # in many papers this value really
computeIntegral = 0            # is 0.1! This is huge!
plateauSlope = 0.8616

activeForceLengthCurveHack_M = createFiberActiveForceLengthCurve(lce0,lce1,lce2,lce3,minActiveForceLengthValueHack,plateauSlope,curviness,computeIntegral)
activeForceLengthParamsHack=["Approx Active Force Length Curve","$f_{act}(l_{MN})$","$l_{MN}$","$f_{act}$","darkorange","orange","gold"]

#-------------------
#Fiber Passive Force

eZero = 0 
eIso  = 0.7
kLow  = 0.2 
kIso  = 2/(eIso-eZero)
curviness = 0.75 
computeIntegral = 1

l0_Mp=1.0
lf_Mp=1.65

l_MpG=np.linspace(l0_Mp,lf_Mp,npts)
FP=[f_p(x) for x in l_MpG]
passiveForceLengthCurve_DG=[l_MpG,FP,l0_Mp,lf_Mp,f_p(l0_Mp),f_p(lf_Mp)]

passiveForceLengthCurve_M = createFiberForceLengthCurve(eZero,eIso,kLow, kIso,curviness,computeIntegral)
passiveForceLengthParams=["Passive Force Length Curve","$f_{pas}(l_{MN})$","$l_{MN}$","$f_{pas}$","dodgerblue","deepskyblue","c","midnightblue"]

#--------------------------
#Fiber Force Velocity Curve
                              
fmaxE=1.4 
dydxC=0 
dydxNearC=0.15
dydxIso=5.0 
dydxE=0.1 
dydxNearE=0.15
flag_smoothenNonZeroDyDxC=0 
flag_usingOctave=0
fvAtHalfVMax=0.15 
concCurviness=0.7 
eccCurviness=0.9 
computeIntegral=0

v0_M=-1
vf_M=1

v_MG=np.linspace(v0_M,vf_M,npts)
FV=[f_v(x) for x in v_MG]
fiberForceVelocityCurve_DG=[v_MG,FV,v0_M,vf_M,min(FV),max(FV)]

fiberForceVelocityCurve_M=createFiberForceVelocityCurve2018(fmaxE,dydxE,dydxC,flag_smoothenNonZeroDyDxC,dydxNearE,fvAtHalfVMax,eccCurviness)
fiberForceVelocityCurveHack_M = createFiberForceVelocityCurve2018(fmaxE,dydxNearE,dydxNearC,flag_smoothenNonZeroDyDxC,dydxNearE,fvAtHalfVMax,eccCurviness)
fiberForceVelocityInverseCurveHack_M = createInverseBezierCurve(fiberForceVelocityCurveHack_M)
fiberForceVelocityParams=["Force Velocity Curve","$f_v(v_{MN})$","$v_{MN}$","$f_v$","darkgreen","green","limegreen"]
fiberForceVelocityParamsHack=["Approx Force Velocity Curve","$f_v(v_{MN})$","$v_{MN}$","$f_v$","darkolivegreen","olive","yellowgreen"]
fiberForceVelocityInverseParamsHack=["Approx Force Velocity Inverse Curve","$f^{inv}_v(v_{MN})$","$v_{MN}$","$f^{inv}_v$","seagreen","mediumseagreen","mediumaquamarine"]

#-------------------------
#Tendon Force Length Curve
kT=35
c1=0.200
c2=0.995
c3=0.250
eIso= 0.049
kIso = 1.375/eIso 
fToe = 2.0/3.0
curviness= 0.5 
computeIntegral = 1

l0_T=(1/kT)*log(c3/c1)+c2
lf_T=1.05

l_TG=np.linspace(l0_T,lf_T,npts)
FT=[f_t(x) for x in l_TG]
tendonForceLengthCurve_DG=[l_TG,FT,l0_T,lf_T,min(FT),max(FT)]

tendonForceLengthCurve_M = createTendonForceLengthCurve(eIso, kIso,fToe, curviness,computeIntegral)
tendonForceLengthParams=["Tendon Force Length Curve","$f_T(l_{TN})$","$l_{TN}$","$f_T$","purple","darkviolet","violet","indigo"]

curveParamVector=[activeForceLengthCurve_M,passiveForceLengthCurve_M,fiberForceVelocityCurve_M,tendonForceLengthCurve_M,activeForceLengthCurveHack_M,fiberForceVelocityCurveHack_M,fiberForceVelocityInverseCurveHack_M]
curveSampleParams=[activeForceLengthParams,passiveForceLengthParams,fiberForceVelocityParams,tendonForceLengthParams,activeForceLengthParamsHack,fiberForceVelocityParamsHack,fiberForceVelocityInverseParamsHack]
DeGrooteModel=[activeForceLengthCurve_DG,passiveForceLengthCurve_DG,fiberForceVelocityCurve_DG,tendonForceLengthCurve_DG]

#Model comparison
ModelCurvesComparison(DeGrooteModel,curveParamVector,curveSampleParams,npts,save)







