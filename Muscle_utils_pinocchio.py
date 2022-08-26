import numpy as np
import sympy
from sympy import symbols, Eq, solve
from DeGroote_Muscle_Utils import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


t_simulation = 0.01

l0_M = 0.55
ls_T = 0.55
alpha0 = 0.209
activation = 0.8
v_max = 10

# Solving for the muscle length using differential equation then calculating FM
def l_MN_calculation(l_MT, l_MN):
    sol = solve_ivp(lambda t, l_MN: ode_DG(t, l_MN, l_MT, activation, l0_M, alpha0, ls_T, v_max), [0, t_simulation],[l_MN], methode="Radau", rtol=1e-6, atol=1e-6,t_eval=np.linspace(0, t_simulation, int(t_simulation * 200)), dense_output=True)
    l_MN = sol.y[0][-1]
    FM = FM_calc(l_MN, activation)
    return FM

def pos_inse(position_second_joint, distance_insertion_to_second_joint, position_end_effector, length_2nd_joint_to_end_effector):
    x_i = position_second_joint[0] + distance_insertion_to_second_joint*(position_end_effector[0]-position_second_joint[0])/length_2nd_joint_to_end_effector
    y_i = position_second_joint[1] + distance_insertion_to_second_joint*(position_end_effector[1]-position_second_joint[1])/length_2nd_joint_to_end_effector
    z_i = position_second_joint[2] + distance_insertion_to_second_joint*(position_end_effector[2]-position_second_joint[2])/length_2nd_joint_to_end_effector
    position_insertion = np.array([x_i, y_i, z_i])
    return position_insertion