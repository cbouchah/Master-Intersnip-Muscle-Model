from __future__ import print_function
from pinocchio.visualize import MeshcatVisualizer
import numpy as np
from numpy.linalg import norm, pinv
from os.path import dirname, join, abspath
import pinocchio as pin
import time
import hppfcl

def display_sphere(des_SE3):
    model2 = pin.Model()
    geom_model = pin.GeometryModel()
    geometries = [hppfcl.Sphere(0.1)]
    for i, geom in enumerate(geometries):
        geom_obj = pin.GeometryObject("obj{}".format(i), 0, 0, geom, des_SE3)
        color = np.random.uniform(0, 1, 4)
        color[3] = 1
        geom_obj.meshColor = color
        geom_model.addGeometryObject(geom_obj)
    viz2 = MeshcatVisualizer(model2, geom_model, geom_model)
    viz2.initViewer()
    viz2.loadViewerModel()
    viz2.display(np.zeros(0))
    return viz2

des_tool_place = pin.SE3()
des_tool_place.translation = np.array([1, 1, 1])
des_tool_place.rotation= pin.utils.rpyToMatrix(np.pi, np.pi/4, np.pi)
des_root_place = pin.SE3(np.eye(3),np.array([0,0,0]))

viz2 = display_sphere(des_tool_place)

pinocchio_model_dir = join("/home/cbouchah/devel", "pinocchio/models")
model_path = join(pinocchio_model_dir, "example-robot-data/robots")
mesh_dir = pinocchio_model_dir
urdf_filename = "romeo_small.urdf"
urdf_model_path = join(join(model_path, "romeo_description/urdf"), urdf_filename)
model, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_model_path, mesh_dir, pin.JointModelFreeFlyer())
data, collision_data, visual_data = pin.createDatas(model, collision_model, visual_model)
viz = MeshcatVisualizer(model, collision_model, visual_model)
viz.initViewer(viz2.viewer)
viz.loadViewerModel(color = [1, 2, 1.5, 2.1])
time.sleep(1)

JOINT_ID = 31
IDX_BASIS = 1
q= pin.neutral(model)
eps    = 1e-4
DT     = 1e-2
viz.display(q)
time.sleep(5)

i=0
j=0

tool_frame_id = model.getFrameId("r_gripper")
root_frame_id = model.getFrameId("l_gripper")

while True :
    while True:
        pin.forwardKinematics(model, data, q)
        pin.updateFramePlacements(model, data)
        current_tool = data.oMf[tool_frame_id]
        err1 = pin.log(current_tool.actInv(des_tool_place)).vector
        if norm(err1) < eps:
            success = True
            break
        if i >= 400:
            success = False
            break
        J1 = pin.computeFrameJacobian(model, data, q, tool_frame_id)
        vq1 = pinv(J1).dot(err1)
        i += 1
        q = pin.integrate(model, q, vq1 * DT)
        print(q)
        viz.display(q)
        if not i % 10:
            print('%d: error1 = %s' % (i, err1.T))


    pin.forwardKinematics(model,data,q)
    pin.updateFramePlacements(model, data)
    current_root = data.oMf[root_frame_id]
    err2 = pin.log(current_root.actInv(des_root_place)).vector
    J1 = pin.computeFrameJacobian(model, data, q, tool_frame_id)
    P1 = np.identity(37) - pinv(J1).dot(J1)

    if norm(err2) < eps:
        success = True
        break
    if j >= 5000:
        success = False
        break
    J2 = pin.computeFrameJacobian(model, data, q, root_frame_id)
    #  vq2 = vq1 + pinv(J2.dot(P1)).dot(err2- J2.dot(vq1))
    vq2 = P1.dot(pinv(J2).dot(err2))
    q = pin.integrate(model,q,vq2*DT)
    print(q)
    viz.display(q)
    j += 1
    if not j % 10:
     print('%d: error2 = %s' % (j, err2.T))



