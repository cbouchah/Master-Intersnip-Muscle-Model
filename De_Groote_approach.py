import time
import sympy
from sympy import symbols, Eq, solve
from DeGroote_Muscle_Utils import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Dimensions of the radius bone and the distance between the elbow joint and of the insertion of the Biceps muscle (Long Head).
l_elbow_to_insertion = 0.055
l_radius = 0.21732

# Coordinates of the origin of the biceps muscle and of the elbow joint.
x_elbow = 1
y_elbow = 1
x_origin = x_elbow
y_origin = y_elbow + l_radius

# Initialize of the muscle and bone representations.
fig_presentation, [[ax_presentation, ax_t_vs_l_MN], [ax_t_differential_vs_l_MN,ax_t_vs_FM]] = plt.subplots(2,2,figsize=(12,12))
lignesb, = ax_presentation.plot([], [], marker='o', c='g', label= "Bones Representation", lw=2)
lignesm, = ax_presentation.plot([], [], marker='o', c='r', label= "Muscle Representation", lw=2)
ax_presentation.set_xlim(0.5, 1.5)
ax_presentation.set_ylim(0.5, 1.5)
lignesm.set_data([x_origin, x_elbow], [y_origin, y_elbow + l_elbow_to_insertion])
lignesb.set_data([x_origin, x_elbow, x_elbow], [y_origin, y_elbow, y_elbow + l_radius])
ax_presentation.set_title('One Joint, One Muscle, Two Bones Representation(2D)')

# Initialize the figure.
l_MN_limitation = (0.5, 1.5)
t_simulation = 0.01
ax_t_differential_vs_l_MN.set_title('Muscle-Tendon Length curve')
ax_t_differential_vs_l_MN.set_xlabel('Time(s)')
ax_t_differential_vs_l_MN.set_ylabel('l_MT')

ax_t_vs_l_MN.set_title('Normalized Muscle Length Curve')
ax_t_vs_l_MN.set_xlabel('Time(s)')
ax_t_vs_l_MN.set_ylabel('Normalised Muscle Length(no unit)')
ax_t_vs_FM.set_title('Muscle Force Curve')
ax_t_vs_FM.set_xlabel('Time(s)')
ax_t_vs_FM.set_ylabel('FM')

# Calculation the coordinates of the insertion and of the wrist then plot it.
def solve_equ_2D_plot(x_elbow, y_elbow, x_origin, y_origin, l_elbow_to_insertion, l_radius, l_MT):
    xc, xb, yc, yb = symbols('xc xb yc yb')
    eq1 = Eq(x_elbow - xc + l_elbow_to_insertion * (xb - x_elbow) / l_radius, 0)
    eq2 = Eq(y_elbow - yc + l_elbow_to_insertion * (yb - y_elbow) / l_radius, 0)
    eq3 = Eq(sympy.sqrt((xb - x_elbow) ** 2 + (yb - x_elbow) ** 2) - l_radius, 0)
    eq4 = Eq(sympy.sqrt((xc - x_origin) ** 2 + (yc - y_origin) ** 2) - l_MT, 0)

    sol = solve((eq1, eq2, eq3, eq4), (xc, xb, yc, yb))
    xc = sol[1][0]
    xb = sol[1][1]
    yc = sol[1][2]
    yb = sol[1][3]
    lignesb.set_data([x_origin, x_elbow, xb], [y_origin, y_elbow, yb])
    lignesm.set_data([x_origin, xc], [y_origin, yc])

# Inputs of the Muscle length calculation.
l0_M = 0.02
c2 = l0_M
ls_T = 0.2
alpha0 = np.pi/6
activation = 0.4
v_max = 10

# Solving for the muscle length using differential equation
def l_MN_calculation(l_MT, l_MN):
    sol = solve_ivp(lambda t,l_MN:ode_DG(t, l_MN, l_MT, activation, l0_M, alpha0, ls_T, v_max),[0, t_simulation],[l_MN],methode="Radau",rtol=1e-6, atol=1e-6, t_eval= np.linspace(0, t_simulation, int(t_simulation*2000)), dense_output=True)
    l_MN = sol.y[0][-1]
    return l_MN

# Turn on interactive plotting and show plot.
plt.ion()
plt.show()

#the first while loop is just for calculating an intial muscle lentgh
t = 0
while True:
    l_MT = calc_l_MT(t,l0_M,ls_T,alpha0,c2)
    # Intial values for l_MT and l_MN
    sol = solve_ivp(lambda t,l_MN:ode_DG(t, l_MN, l_MT, activation, l0_M, alpha0, ls_T, v_max),[0, 0.02],[1],methode="Radau",rtol=1e-6, atol=1e-6, t_eval= np.linspace(0, 0.02, int(0.2*2000)), dense_output=True)
    l_MN = sol.y[0][-1]
    FM = FM_calc(l_MN,activation)
    t+=0.01
    if t>0.05:
        success = False
    break

#while loop to calculate at each time all the parameters
while True:
    if t>0.7:
        success = False
        print("The figure will disappear in a minute and a half.")
        time.sleep(900)
        break
    l_MN_old = l_MN
    FM_old = FM
    l_MT_old = l_MT
    solve_equ_2D_plot(x_elbow, y_elbow, x_origin, y_origin, l_elbow_to_insertion, l_radius, l_MT)
    l_MN = l_MN_calculation(l_MT, l_MN)
    FM = FM_calc(l_MN,activation)
    l_MT = calc_l_MT(t,l0_M,ls_T,alpha0,c2)
    ax_t_differential_vs_l_MN.plot([t,t+0.01], [l_MT_old, l_MT], color='blue')
    ax_t_vs_l_MN.plot([t,t+0.01], [l_MN_old, l_MN], color='blue')
    ax_t_vs_FM.plot([t,t+0.01], [FM_old, FM], color='blue')
    t+=0.01
    fig_presentation.canvas.get_tk_widget().update()