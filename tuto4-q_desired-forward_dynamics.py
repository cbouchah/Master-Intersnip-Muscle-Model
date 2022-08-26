from __future__ import print_function
import example_robot_data as robex
import pinocchio as pin
import numpy as np
import time

name = "talos_arm"
robot = robex.load(name, rootNodeName="talos_arm")

Viewer = pin.visualize.MeshcatVisualizer
viz = Viewer(robot.model, robot.collision_model, robot.visual_model)
viz.initViewer(loadModel=True)
q = np.array([0, 0, 0.2, 0, 0, 0.3,0])*np.pi
qdes = np.array([0, 0, 0.5, 0, 0, 0.8,0])*np.pi
v = np.zeros(robot.model.nv)
a = np.zeros(robot.model.nv)
dt = 1e-4
torq = np.zeros(robot.model.nv)
Kp = 110
Kv = 2 * np.sqrt(Kp)
viz.display(qdes)
time.sleep(4)

while True :
    torq = -Kp * (q - qdes) - Kv * v
    # b = pin.nle(robot.model, robot.data, q, v)
    # M = pin.crba(robot.model, robot.data, q)
    # a_free = pinv(M) .dot(torq - b)
    a_free = pin.aba(robot.model, robot.data, q, v, torq)
    v += a_free * dt
    q = pin.integrate(robot.model, q, v * dt)
    viz.display(q)
