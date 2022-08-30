# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/createDefaultNormalizedMuscleCurves.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine
from Mathematical_expressions_for_muscle_tendon_characteristics_De_Groote import F_MN
from calcBezierYFcnXCurveSampleVector import calcBezierYFcnXCurveSampleVector
import matplotlib.pyplot as plt

#Model comparison
def ModelCurvesComparison(DeGrooteModel,curveParamVector,curveSampleParams,npts,save):
    
    """Fonction qui va tracer les différentes courbes des deux modèles pour pouvoir les comparer
    
    -DeGrooteModel liste qui regroupe l'ensemble des forces que l'on veut tracer
    et chaque élements DeGrooteModel[i]=(x,f(x),xmin,xmax,ymin,ymax)
    
    -curveParamVector liste qui regroupe l'ensemble des forces que l'on veut tracer
    et chaque élements curveParamVector[i]=(x,y,dydx,d2ydx2,d3ydx3,intYdx)
    
    -curveSampleParams liste qui regroupe l'ensemble des paramètres/légendes que l'on va 
    utiliser pour tracer les courbes, 
    curveSampleParams[i]=(title,labelfonc,xlabel,ylabel,c1,c2,c3, (c4) )
    
    -npts int, le nb de point que l'on choisit
    
    -save bol, True si l'on souhaite sauvegarder les courbes, False sinon"""
    
    for i in range(len(DeGrooteModel)):
        curveSample_M = calcBezierYFcnXCurveSampleVector(curveParamVector[i],npts)
        
        xmin =min(min(curveSample_M[0]),DeGrooteModel[i][2])
        xmax = max(max(curveSample_M[0]),DeGrooteModel[i][3])
        ymin = min(min(curveSample_M[1]),DeGrooteModel[i][4])
        ymax = max(max(curveSample_M[1]),DeGrooteModel[i][5])
        yDelta = ymax-ymin
        
        xV   = curveSample_M[0]
        yV   = curveSample_M[1]
        
        plt.figure(curveSampleParams[i][0]+" - Comparison")
        plt.plot(DeGrooteModel[i][0],DeGrooteModel[i][1],label=curveSampleParams[i][1]+" De Groote Model",color=curveSampleParams[i][6])
        plt.plot(xV,yV,label=curveSampleParams[i][1]+" Millard Model",color=curveSampleParams[i][4])
        plt.xlabel(curveSampleParams[i][2])
        plt.ylabel(curveSampleParams[i][3])
        plt.title(curveSampleParams[i][0]+" - Comparison")
        plt.legend()
        plt.axis ([xmin,xmax,ymin-0.1*yDelta,ymax+0.1*yDelta])
        plt.axhline(color="black")
        plt.axvline(color="black")
        if save==True:
            plt.savefig(curveSampleParams[i][0]+" - Comparison")
        plt.show()
        
    """ Pour une meilleure visualisation de la force musculaire, 
    on représente ici la force active et passive sur un même graphique:"""
        
    curveSample1 = calcBezierYFcnXCurveSampleVector(curveParamVector[0],npts)
    curveSample2 = calcBezierYFcnXCurveSampleVector(curveParamVector[1],npts)
    
    xmin = min(min(curveSample1[0]),min(curveSample2[0]),DeGrooteModel[0][2],DeGrooteModel[1][2])
    xmax = max(max(curveSample1[0]),max(curveSample2[0]),DeGrooteModel[0][3],DeGrooteModel[1][3])
    ymin = min(min(curveSample1[1]),min(curveSample2[1]),DeGrooteModel[0][4],DeGrooteModel[1][4])
    ymax = max(max(curveSample1[1]),max(curveSample2[1]),DeGrooteModel[0][5],DeGrooteModel[1][5])
    yDelta = ymax-ymin
    
    xV1   = curveSample1[0]
    yV1   = curveSample1[1]
    xV2   = curveSample2[0]
    yV2   = curveSample2[1]
    
    plt.figure("Active and Passive Force Length Curve - Comparison model")
    plt.plot(xV1,yV1,label="$f_{act}(l_{MN})$ Millard Model",color="r")
    plt.plot(DeGrooteModel[0][0],DeGrooteModel[0][1],label=curveSampleParams[0][1]+" De Groote Model",color=curveSampleParams[0][6])
    plt.plot(xV2,yV2,label="$f_{pas}(l_{MN})$ Millard Model",color="deepskyblue")
    plt.plot(DeGrooteModel[1][0],DeGrooteModel[1][1],label=curveSampleParams[1][1]+" De Groote Model",color=curveSampleParams[1][6])
    plt.xlabel("$l_{MN}$")
    plt.ylabel("$f_{act}$ and $f_{pas}$")
    plt.title("Active and Passive Force Length Curve - Comparison model")
    plt.legend()
    plt.axis ([xmin,xmax,ymin-0.1*yDelta,ymax+0.1*yDelta])
    plt.axhline(color="black")
    plt.axvline(color="black")
    if save==True:
        plt.savefig("Active and Passive Force Length Curve - Comparison model")
    plt.show()
        
    return None

