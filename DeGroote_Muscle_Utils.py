from math import pi, sqrt, sin, sinh
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
    return (np.exp((kpe*(l_MN-1))/e0)-1)/(np.exp(kpe)-1)

#The velocity force-length characteristic : f__v
def f_v(v_MN):
    v = d1*np.log((d2*v_MN+d3)+np.sqrt(((d2*v_MN+d3)**2)+1))+d4
    return v

#muscle-tendon length variation through time
def calc_l_MT(t,l0_M,ls_T,alpha0,c2):
    c1 = ls_T + l0_M*(np.cos(alpha0))
    return c1 + c2*(np.sin(2*np.pi*t))

#muscle force
def FM_calc(l_MN, a):
    FM = a*f_a(l_MN)+f_p(l_MN)
    return FM

#differential equation inspired by DeGroote
def ode_DG(t, l_MN, l_MT, a,l0_M,alpha0,ls_T,vmax_M):

    assert 0<alpha0 and alpha0<=pi/2
    if l_MN**2-sin(alpha0)**2 <0:
        sqt = sqrt(sin(alpha0)**2-l_MN**2)
    else:
        sqt = sqrt(l_MN**2-sin(alpha0)**2)

    l_T= l_MT -l0_M*(np.absolute(sqt))
    CosA=sqt/l_MN

    A=a*f_a(l_MN)
    B=f_t(l_T/ls_T)/CosA
    C=f_p(l_MN)
    FM=(B-C)/A
    if FM>d4-710.47*d1:
        FM=d4-710.47*d1
    if FM<d4+710.47*d1:
        FM=d4+710.47*d1
    r=(1/d1)*(FM-d4)
    finv_v=(1/d2)*sinh(r)-d3
    return  (vmax_M/l0_M)*finv_v
