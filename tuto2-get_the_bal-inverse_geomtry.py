import pinocchio as pin
import numpy as np
from scipy.optimize import fmin_bfgs, fmin_slsqp
from os.path import dirname, join, abspath
from pinocchio.visualize import MeshcatVisualizer
import time
import hppfcl

#steps to load romeo robot
pinocchio_model_dir = join("/home/cbouchah/devel", "pinocchio/models")
model_path = join(pinocchio_model_dir, "example-robot-data/robots")
mesh_dir = pinocchio_model_dir
urdf_filename = "romeo_small.urdf"
urdf_model_path = join(join(model_path, "romeo_description/urdf"), urdf_filename)
model, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_model_path, mesh_dir, pin.JointModelFreeFlyer())
data, collision_data, visual_data = pin.createDatas(model, collision_model, visual_model)

#steps to load a sphere
model2 = pin.Model()
geom_model = pin.GeometryModel()
geometries = [hppfcl.Sphere(0.1)]
for i, geom in enumerate(geometries):
    placement = pin.SE3(np.eye(3), np.array([0.33, -0.25, 0.6]))
    geom_obj = pin.GeometryObject("obj{}".format(i), 0, 0, geom, placement)
    color = np.random.uniform(0, 1, 4)
    color[3] = 1
    geom_obj.meshColor = color
    geom_model.addGeometryObject(geom_obj)

viz2 = MeshcatVisualizer(model2, geom_model, geom_model)
viz2.initViewer()
viz2.loadViewerModel()
viz2.display(np.zeros(0))

#steps to load the robot and the sphere
viz = MeshcatVisualizer(model, collision_model, visual_model)
viz.initViewer(viz2.viewer)
viz.loadViewerModel(color = [1, 2, 1.5, 2.1])

#cost function inputs (q and desire position) output (distance between the q given and the desire)
def cost(q, des1, des2, des3):
    qtot=np.concatenate((np.array([0, 0, 0.840252, 0, 0, 0, 1]),q),axis=None)
    pin.forwardKinematics(model, data, qtot)
    pos = data.oMi[31].translation
    return np.linalg.norm(pos-np.array([des1, des2, des3]))

#function to display q
def disp(x):
    time.sleep(1)
    xtot=np.concatenate((np.array([0, 0, 0.840252, 0, 0, 0, 1]),x),axis=None)
    viz.display(xtot)
    print(xtot)

x0 = np.zeros(31)
disp(x0)

# to minimize the distance between the end effector and the desire position(sphere)
xopt_bfgs = fmin_bfgs(cost, x0, args=np.array([0.33, -0.25, 0.6]),callback=disp)


