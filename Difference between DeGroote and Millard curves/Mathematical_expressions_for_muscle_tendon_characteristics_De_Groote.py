# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------------------------------------------------
#DE GROOTE MODEL
#From "Evaluation of Direct Collocation Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem"
#Written by: Friedl De Groote, Allison L. Kinney, Anil V. Rao, Benjamin J. Fregly
#--> All numerical values ​​are available on the Online Supplement
#----------------------------------------------------------------------------------------------------------------------
#Edited by christian bou chahine

import numpy as np

#Parameters od the Hill model characteristics 
#Tendon force-lengh
kT=35
c1=0.200
c2=0.995
c3=0.250
#-------------------------
#Active muscle force-length
b11=0.815
b21=1.055
b31=0.162
b41=0.063
b12=0.433
b22=0.717
b32=-0.030
b42=0.200
b13=0.100
b23=1.000
b33=0.354
b43=0.000
#-------------------------
#Passive muscle force-length
kpe=4.0
e0=0.6
#-------------------------
#Muscle force-velocity
d1=-0.318
d2=-8.149
d3=-0.374
d4=0.886
#-------------------------

#The tendon force-length characteristic : f_t
#f_t=lambda l_T: c1*np.exp(kT*(l_T-c2))-c3 
def f_t(l_TN):
    if c1*np.exp(kT*(l_TN-c2))-c3<0 :
        return 0
    else:
        return c1*np.exp(kT*(l_TN-c2))-c3 
    
#The active force-length characteristic : f_a
def f_a(l_MN): #with l_M normalized fiber length
    S1=b11*np.exp((-0.5*(l_MN-b21)**2)/(b31+b41*l_MN)**2)
    S2=b12*np.exp((-0.5*(l_MN-b22)**2)/(b32+b42*l_MN)**2)
    S3=b13*np.exp((-0.5*(l_MN-b23)**2)/(b33+b43*l_MN)**2)
    return S1+S2+S3 

#The passive force-length characteristic : f__p
def f_p(l_MN):
    if l_MN<1:
        return 0
    else:
        return (np.exp((kpe*(l_MN-1))/e0)-1)/(np.exp(kpe)-1)

#The total of force-length characteristic : f_a+f_p
F_MN=lambda l_MN : f_a(l_MN)+f_p(l_MN)


#The force-velocity characteristic : f_v
f_v=lambda v_MN: d1*np.log((d2*v_MN+d3)+np.sqrt(((d2*v_MN+d3)**2)+1))+d4

    


