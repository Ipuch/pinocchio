"""Demonstrate the spline-joint regression algorithm ``pin.fitSplineJoint``.

Companion of ``spline-joint.py``: there the spline control frames are written by
hand; here we *recover* them. We sample SE(3) placements along a target
trajectory (a parabola ``z = y^2`` in the y-z plane) and fit a spline joint to
them by Gauss-Newton regression. The fitted joint reproduces the parabola, so a
box dropped on it slides under gravity exactly as in ``spline-joint.py``.

The Meshcat scene shows the regression input (green data placements), the
recovered control frames (blue dots) and the box sliding on the fitted spline.
"""

import time

import meshcat.geometry as mg
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer


class PointTracer:
    """Trace a fixed point in Meshcat, keeping the dotted trail."""

    def __init__(self, viz):
        self._trail = []
        self._node = viz.viewer["solid_pose"]

    def add(self, solid_pose):
        self._trail.append(solid_pose.translation.copy())
        pts = np.asarray(self._trail, dtype=np.float32).T
        self._node.set_object(
            mg.Line(mg.PointsGeometry(pts), mg.LineBasicMaterial(color=0xFF3030))
        )


def draw_points(viz, name, points, color, size):
    """Draw a set of points as a Meshcat point cloud."""
    pts = np.asarray(points, dtype=np.float32).T
    viz.viewer[name].set_object(
        mg.Points(mg.PointsGeometry(pts), mg.PointsMaterial(size=size, color=color))
    )


def parabola_data(n_samples=60, y_min=-1.0, y_max=1.0):
    """Sample SE(3) placements along the parabola z = y^2 in the y-z plane.

    Returns (data_frames, q_data) with q_data = y, so the recovered joint is
    driven directly in y over [y_min, y_max].
    """
    ys = np.linspace(y_min, y_max, n_samples)
    data_frames = [pin.SE3(np.eye(3), np.array([0.0, y, y * y + 0.1* np.exp(y) + 0.1 * np.random.rand()])) for y in ys]
    return data_frames, ys


def sim_loop(model, dt=0.01, nsteps=800):
    qs = [np.array([1.0])]
    vs = [np.array([0.0])]
    data = model.createData()
    for i in range(nsteps):
        q = qs[i]
        v = vs[i]
        tau = -1.0 * v  # a little bit of damping
        a = pin.aba(model, data, q, v, tau)
        vnext = v + dt * a
        qnext = pin.integrate(model, q, dt * vnext)
        qs.append(qnext)
        vs.append(vnext)
    return qs, vs


def main():
    # 1. Regression input: placements sampled along the target parabola.
    data_frames, q_data = parabola_data()

    # 2. Recover the spline joint by Gauss-Newton regression on those placements.
    settings = pin.SplineRegressionSettings()
    spline_joint = pin.fitSplineJoint(data_frames, q_data, 10, 8, settings)
    print(
        f"regression: {settings.iter} iterations, "
        f"residual {settings.absolute_residual:.2e}"
    )
    print(
        f"recovered {len(spline_joint.ctrlFrames)} control frames, "
        f"q in [{spline_joint.min_q:.1f}, {spline_joint.max_q:.1f}]"
    )

    # 3. Same setup as spline-joint.py: drop a box on the fitted parabola.
    model = pin.Model()
    joint_id = model.addJoint(0, spline_joint, pin.SE3.Identity(), "spline-joint")
    box_inertia = pin.Inertia.FromBox(1.0, 1.0, 1.0, 1.0)
    model.appendBodyToJoint(joint_id, box_inertia, pin.SE3.Identity())

    try:
        viz = MeshcatVisualizer(model, None, pin.GeometryModel())
        viz.initViewer(open=True)
        viz.loadViewerModel()
    except ImportError as e:
        print("Error while initializing the viewer.")
        print(e)
        return

    # Show the regression input and the recovered control frames.
    draw_points(
        viz, "data_frames", [X.translation for X in data_frames], 0x30C030, 0.02
    )
    draw_points(
        viz,
        "control_frames",
        [C.translation for C in spline_joint.ctrlFrames],
        0x3060FF,
        0.04,
    )

    solid_frame = viz.viewer["solid_frame"]
    solid_frame.set_object(mg.triad(0.2))
    tracer = PointTracer(viz)

    dt = 0.01
    qs, _ = sim_loop(model, dt)
    for q in qs:
        viz.display(q)
        solid_pose = viz.data.oMi[joint_id]
        tracer.add(solid_pose)
        solid_frame.set_transform(solid_pose.homogeneous)
        time.sleep(dt)


if __name__ == "__main__":
    main()
