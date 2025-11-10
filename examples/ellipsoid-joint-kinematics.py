"""
Ellipsoid Joint Kinematics Example

This example demonstrates the use of the ellipsoid joint in Pinocchio.
The ellipsoid joint constrains a point to move on the surface of an ellipsoid
while allowing full 3D rotational freedom (SO(3)).

The joint has 3 DOF corresponding to 3 rotation angles that determine both:
- The position on the ellipsoid surface
- The orientation of the body

This is useful for modeling:
- Joints with coupled translation and rotation
- Contact points on ellipsoidal surfaces
- Biomechanical joints with ellipsoidal constraints
"""

import pinocchio as pin
import numpy as np

# For visualization (optional)
try:
    from pinocchio.visualize import MeshcatVisualizer
    VISUALIZE = True
except ImportError:
    print("MeshcatVisualizer not available. Skipping visualization.")
    VISUALIZE = False


def create_ellipsoid_robot(radius_a=0.5, radius_b=0.3, radius_c=0.2):
    """
    Create a simple robot with an ellipsoid joint.
    
    Parameters
    ----------
    radius_a : float
        Ellipsoid radius along x-axis
    radius_b : float
        Ellipsoid radius along y-axis  
    radius_c : float
        Ellipsoid radius along z-axis
        
    Returns
    -------
    model : pin.Model
        The robot model
    """
    model = pin.Model()
    
    # Create ellipsoid joint
    joint_name = "ellipsoid"
    joint_id = model.addJoint(
        0,  # parent joint id (0 = universe)
        pin.JointModelEllipsoid(radius_a, radius_b, radius_c),
        pin.SE3.Identity(),
        joint_name
    )
    
    # Add a body with some inertia
    mass = 1.0
    lever = np.array([0.0, 0.0, 0.1])
    inertia = pin.Inertia(
        mass,
        lever,
        np.diag([0.01, 0.01, 0.01])
    )
    
    body_placement = pin.SE3.Identity()
    body_placement.translation = np.array([0.0, 0.0, 0.2])
    
    model.appendBodyToJoint(joint_id, inertia, body_placement)
    
    return model


def compute_kinematics_example():
    """Demonstrate forward kinematics with the ellipsoid joint."""
    print("=" * 60)
    print("ELLIPSOID JOINT KINEMATICS EXAMPLE")
    print("=" * 60)
    
    # Create model with ellipsoid radii
    radius_a, radius_b, radius_c = 0.5, 0.3, 0.2
    model = create_ellipsoid_robot(radius_a, radius_b, radius_c)
    data = model.createData()
    
    print(f"\nEllipsoid parameters:")
    print(f"  radius_a (x-axis): {radius_a}")
    print(f"  radius_b (y-axis): {radius_b}")
    print(f"  radius_c (z-axis): {radius_c}")
    
    # Test different configurations
    configs = [
        ("Zero configuration", np.array([0.0, 0.0, 0.0])),
        ("Rotation around x", np.array([np.pi/4, 0.0, 0.0])),
        ("Rotation around y", np.array([0.0, np.pi/4, 0.0])),
        ("Rotation around z", np.array([0.0, 0.0, np.pi/4])),
        ("Combined rotations", np.array([np.pi/6, np.pi/4, np.pi/3])),
    ]
    
    for name, q in configs:
        pin.forwardKinematics(model, data, q)
        pin.updateFramePlacements(model, data)
        
        placement = data.oMi[1]
        position = placement.translation
        
        print(f"\n{name}:")
        print(f"  q = {q}")
        print(f"  Position on ellipsoid: {position}")
        print(f"  Distance from origin: {np.linalg.norm(position):.4f}")
        
        # Verify point is on ellipsoid surface
        normalized_pos = position / np.array([radius_a, radius_b, radius_c])

        # checking the quadratic form of the ellipsoid equation
        on_surface = np.abs(np.linalg.norm(normalized_pos) - 1.0) < 1e-10
        print(f"  On ellipsoid surface: {on_surface}")


def compute_dynamics_example():
    """Demonstrate dynamics computation with the ellipsoid joint."""
    print("\n" + "=" * 60)
    print("ELLIPSOID JOINT DYNAMICS EXAMPLE")
    print("=" * 60)
    
    # Create model
    model = create_ellipsoid_robot(0.5, 0.3, 0.2)
    data = model.createData()
    
    # Configuration, velocity, and acceleration
    q = np.array([np.pi/6, np.pi/4, np.pi/3])
    v = np.array([0.1, 0.2, 0.3])
    a = np.array([0.01, 0.02, 0.03])
    
    print(f"\nConfiguration: q = {q}")
    print(f"Velocity:      v = {v}")
    print(f"Acceleration:  a = {a}")
    
    # Compute dynamics with RNEA (Recursive Newton-Euler Algorithm)
    # Important Note: the generalized forces dimension matches the number of DOFs (3 for ellipsoid joint)
    #   but they are not pure torques due to the joint's nature and because the last
    #   The last step of rnea is tau = S.T * f
    #   where S is the joint motion subspace, and f the spatial force of the joint
    #   So the resulting tau vector contains generalized forces that
    #   include contributions from both rotational and translational effects.
    tau = pin.rnea(model, data, q, v, a)
    print(f"\nRequired torques (RNEA): τ = {tau}")

    # separate the generalized forces for interpretability from the spatial forces
    S = data.joints[1].S
    f = np.array(data.f[1])
    pure_torques = S[3:, :].T @ f[3:]  # Extract rotational part
    translations = S[:3, :].T @ f[:3]  # Extract translational part
    print(f"  Pure rotational torques: {pure_torques}")
    print(f"  Translational forces: {translations}")

    # add them back to recover the original tau
    tau_reconstructed = np.zeros(3)
    tau_reconstructed += pure_torques
    tau_reconstructed += translations
    print(f"  Reconstructed τ: {tau_reconstructed}")



def visualize_ellipsoid_motion():
    """Visualize the ellipsoid joint motion in MeshCat."""
    if not VISUALIZE:
        print("\nVisualization skipped (MeshcatVisualizer not available)")
        return

    if not hasattr(pin, 'hppfcl'):
        print("\nVisualization skipped (hpp-fcl not available)")
        print("Install with: conda install -c conda-forge hpp-fcl")
        print("Then recompile Pinocchio with -DBUILD_WITH_HPP_FCL_SUPPORT=ON")
        return

    print("\n" + "=" * 60)
    print("ELLIPSOID JOINT VISUALIZATION")
    print("=" * 60)

    # Create model
    radius_a, radius_b, radius_c = 0.5, 0.3, 0.2
    model = create_ellipsoid_robot(radius_a, radius_b, radius_c)

    # Add visual geometry
    try:
        geom_model = pin.GeometryModel()

        # 1. Add the ELLIPSOID SURFACE as a visual object
        ellipsoid_shape = pin.hppfcl.Ellipsoid(radius_a, radius_b, radius_c)
        ellipsoid_geom = pin.GeometryObject(
            "ellipsoid_surface",
            0,  # Universe frame
            ellipsoid_shape,
            pin.SE3.Identity()
        )
        ellipsoid_geom.meshColor = np.array([0.8, 0.8, 0.8, 0.3])  # Semi-transparent gray
        geom_model.addGeometryObject(ellipsoid_geom)

        # 2. Add a small sphere to show the contact point on the ellipsoid
        contact_sphere = pin.hppfcl.Sphere(0.03)
        contact_geom = pin.GeometryObject(
            "contact_point",
            model.getJointId("ellipsoid"),
            contact_sphere,
            pin.SE3.Identity()  # At the joint frame (on the ellipsoid surface)
        )
        contact_geom.meshColor = np.array([1.0, 0.0, 0.0, 1.0])  # Red
        geom_model.addGeometryObject(contact_geom)

        # 3. Add a box to visualize the body orientation
        box = pin.hppfcl.Box(0.1, 0.1, 0.3)
        box_geom = pin.GeometryObject(
            "body_visual",
            model.getJointId("ellipsoid"),
            box,
            pin.SE3(np.eye(3), np.array([0.0, 0.0, 0.15]))
        )
        box_geom.meshColor = np.array([0.2, 0.6, 1.0, 1.0])  # Blue
        geom_model.addGeometryObject(box_geom)

        # Initialize visualizer
        viz = MeshcatVisualizer(model, geom_model, geom_model)
        viz.initViewer(open=True)
        viz.loadViewerModel()

        print("\n✓ Visualization initialized!")
        print("  - Gray ellipsoid: The constraint surface")
        print("  - Red sphere: Contact point on the surface")
        print("  - Blue box: The rigid body attached to the joint")
        print("\nAnimating ellipsoid joint motion...")
        print("Open http://127.0.0.1:7000/static/ in your browser to see the animation.")

        # Animate through different configurations
        import time
        t = 0.0
        dt = 0.02

        for i in range(500):
            # Create a smooth trajectory on the ellipsoid
            q = np.array([
                0.8 * np.sin(1.0 * t),
                0.6 * np.sin(1.2 * t + 1.0),
                0.4 * np.sin(0.8 * t + 0.5)
            ])

            viz.display(q)
            time.sleep(dt)
            t += dt

            if i % 50 == 0:
                print(f"  Frame {i}/500 - Configuration: [{q[0]:.3f}, {q[1]:.3f}, {q[2]:.3f}]")

        # create extra motion that slides along the principal axes q0
        for angle in np.linspace(0, 2 * np.pi, 100):
            q = np.array([
                0.5 * np.cos(angle),
                0.0,
                0.0
            ])
            viz.display(q)
            time.sleep(dt)

        # q1 axis
        for angle in np.linspace(0, 2 * np.pi, 100):
            q = np.array([
                0.0,
                0.3 * np.cos(angle),
                0.0
            ])
            viz.display(q)
            time.sleep(dt)

        # q2 axis
        for angle in np.linspace(0, 2 * np.pi, 100):
            q = np.array([
                0.0,
                0.0,
                0.2 * np.cos(angle)
            ])
            viz.display(q)
            time.sleep(dt)

        print("\n✓ Animation completed!")

    except Exception as e:
        print(f"\nVisualization error: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all examples."""
    compute_kinematics_example()
    compute_dynamics_example()
    visualize_ellipsoid_motion()
    
    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
