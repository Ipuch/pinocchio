//
// Copyright (c) 2025 INRIA
//

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/multibody/joint/joints.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/frames.hpp"

#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace pinocchio;

template<typename D>
void addJointAndBody(
  Model & model,
  const JointModelBase<D> & jmodel,
  const Model::JointIndex parent_id,
  const SE3 & joint_placement,
  const std::string & joint_name,
  const Inertia & Y)
{
  Model::JointIndex idx;

  idx = model.addJoint(parent_id, jmodel, joint_placement, joint_name);
  model.appendBodyToJoint(idx, Y);
}

BOOST_AUTO_TEST_SUITE(JointEllipsoid)

BOOST_AUTO_TEST_CASE(vsSphericalZYX)
{
  using namespace pinocchio;
  typedef SE3::Vector3 Vector3;
  typedef SE3::Matrix3 Matrix3;

  Inertia inertia(1., Vector3(0.5, 0., 0.0), Matrix3::Identity());
  SE3 pos(1);
  pos.translation() = SE3::LinearType(1., 0., 0.);

  // Create both models with the same structure
  Model modelEllipsoid, modelSphericalZYX;
  addJointAndBody(modelEllipsoid, JointModelEllipsoid(0, 0, 0), 0, pos, "ellipsoid", inertia);
  addJointAndBody(modelSphericalZYX, JointModelSphericalZYX(), 0, pos, "spherical", inertia);

  Data dataEllipsoid(modelEllipsoid);
  Data dataSphericalZYX(modelSphericalZYX);

  // Start with ZYX angles
  Eigen::VectorXd q_s(3);
  q_s << 0.5, 1.2, -0.8; // Z=0.5, Y=1.2, X=-0.8
  Eigen::VectorXd qd_s(3);
  qd_s << 0.1, -0.3, 0.7; // ZYX angle velocities
  Eigen::VectorXd qdotdot_s(3);
  qdotdot_s << 0.2, 0.1, -0.1; // ZYX angle accelerations
  // Compute the rotation matrix from ZYX angles
  forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_s, qd_s);
  const Matrix3 & R = dataSphericalZYX.oMi[1].rotation();
  const Motion & spatial_vel_zyx = dataSphericalZYX.v[1];

  // Extract XYZ Euler angles from the rotation matrix
  // For XYZ convention: R = Rx(x) * Ry(y) * Rz(z)
  // We need to solve for x, y, z

  Eigen::Vector3d q_e;

  // From the Ellipsoid rotation matrix structure:
  // R(0,2) = sin(y)
  // R(1,2) = -sin(x)*cos(y)
  // R(2,2) = cos(x)*cos(y)

  double sy = R(0, 2);
  q_e(1) = std::asin(sy); // y angle

  double cy = std::cos(q_e(1));

  if (std::abs(cy) > 1e-6)
  {
    // Not at singularity
    q_e(0) = std::atan2(-R(1, 2) / cy, R(2, 2) / cy); // x angle
    q_e(2) = std::atan2(-R(0, 1) / cy, R(0, 0) / cy); // z angle
  }
  else
  {
    // Gimbal lock - choose x = 0
    q_e(0) = 0.0;
    q_e(2) = std::atan2(R(1, 0), R(1, 1));
  }

  // Get the motion subspace matrices (which give us the Jacobians)
  JointModelSphericalZYX jmodel_s;
  jmodel_s.setIndexes(0, 0, 0);

  JointDataSphericalZYX jdata_s = jmodel_s.createData();
  jmodel_s.calc(jdata_s, q_s);

  JointModelEllipsoid jmodel_e(0, 0, 0);
  jmodel_e.setIndexes(0, 0, 0);

  JointDataEllipsoid jdata_e = jmodel_e.createData();
  jmodel_e.calc(jdata_e, q_e);

  // The motion subspace S gives us: omega = S * v
  Matrix3 S_s = jdata_s.S.matrix().bottomRows<3>(); // Angular part
  Matrix3 S_e = jdata_e.S.matrix().bottomRows<3>(); // Angular part

  Eigen::Vector3d qd_e = S_e.inverse() * S_s * qd_s;

  Eigen::Vector3d w_s = S_s * qd_s;
  Eigen::Vector3d w_e = S_e * qd_e;

  BOOST_CHECK(w_s.isApprox(w_e));
  
  // Compute forward kinematics with the converted configurations
  forwardKinematics(modelEllipsoid, dataEllipsoid, q_e, qd_e);

  // Getting S with q_e from the three calcs
  JointDataEllipsoid jDataEllipsoidFK = jmodel_e.createData();
  JointDataEllipsoid jDataEllipsoidFK2 = jmodel_e.createData();
  JointDataEllipsoid jDataEllipsoidFK3 = jmodel_e.createData();

  jmodel_e.calc(jDataEllipsoidFK, q_e, qd_e);
  jmodel_e.calc(jDataEllipsoidFK2, q_e);
  jmodel_e.calc(jDataEllipsoidFK3, q_e);
  jmodel_e.calc(jDataEllipsoidFK3, Blank(), qd_e);

  BOOST_CHECK(jDataEllipsoidFK.S.matrix().isApprox(jDataEllipsoidFK2.S.matrix()));
  BOOST_CHECK(jDataEllipsoidFK.S.matrix().isApprox(jDataEllipsoidFK3.S.matrix()));

  JointDataSphericalZYX jDataSphereFK = jmodel_s.createData();
  jmodel_s.calc(jDataSphereFK, q_s, qd_s);

  // Joint-frame velocities (S * v)
  Eigen::Matrix<double, 6, 1> joint_vel_e = jDataEllipsoidFK.S.matrix() * qd_e;
  Eigen::Matrix<double, 6, 1> manual_vel = jDataEllipsoidFK.S.matrix() * jDataEllipsoidFK.joint_v;
  Eigen::Matrix<double, 6, 1> joint_vel_s = jDataSphereFK.S.matrix() * qd_s;

  BOOST_CHECK(dataEllipsoid.v[1].toVector().isApprox(dataSphericalZYX.v[1].toVector()));
  BOOST_CHECK(dataEllipsoid.oMi[1].isApprox(dataSphericalZYX.oMi[1]));

  Matrix3 Sdot_e = jDataEllipsoidFK.Sdot.matrix().bottomRows<3>(); // Angular part

  // The acceleration conversion formula: wdot_s = wdot_e
  // S_s * qdotdot_s + c_s.angular() = Sdot_e * qd_e + S_e * qdotdot_e
  // Solving for qdotdot_e:
  // S_e * qdotdot_e =  c_s.angular()+ S_s * qdotdot_s - Sdot_e * qd_e
  Eigen::Vector3d qdotdot_e =
    S_e.inverse() * (S_s * qdotdot_s + jDataSphereFK.c.angular() - Sdot_e * qd_e);

  // Verify angular accelerations match
  Eigen::Vector3d wdot_s = jDataSphereFK.c.angular() + S_s * qdotdot_s;
  Eigen::Vector3d wdot_e = Sdot_e * qd_e + S_e * qdotdot_e;
  BOOST_CHECK(wdot_s.isApprox(wdot_e));

  forwardKinematics(modelEllipsoid, dataEllipsoid, q_e, qd_e, qdotdot_e);
  forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_s, qd_s, qdotdot_s);

  BOOST_CHECK(dataEllipsoid.a[1].toVector().isApprox(dataSphericalZYX.a[1].toVector()));

  // Test RNEA (Recursive Newton-Euler Algorithm) - spatial forces should match
  rnea(modelEllipsoid, dataEllipsoid, q_e, qd_e, qdotdot_e);
  rnea(modelSphericalZYX, dataSphericalZYX, q_s, qd_s, qdotdot_s);

  BOOST_CHECK(dataEllipsoid.f[1].isApprox(dataSphericalZYX.f[1]));

  // Test ABA (Articulated-Body Algorithm)
  Eigen::VectorXd tau = Eigen::VectorXd::Ones(modelEllipsoid.nv);
  Eigen::VectorXd aAbaEllipsoid =
    aba(modelEllipsoid, dataEllipsoid, q_e, qd_e, dataEllipsoid.tau, Convention::WORLD);

  BOOST_CHECK(dataEllipsoid.ddq.isApprox(qdotdot_e));
}

BOOST_AUTO_TEST_CASE(vsCompositeTxTyTzRxRyRz)
{
  using namespace pinocchio;
  typedef SE3::Vector3 Vector3;
  typedef SE3::Matrix3 Matrix3;

  // Ellipsoid parameters
  double radius_a = 2.0;
  double radius_b = 1.5;
  double radius_c = 1.0;

  Inertia inertia = Inertia::Identity();
  SE3 pos = SE3::Identity();

  // Create Ellipsoid model
  Model modelEllipsoid;
  JointModelEllipsoid jointModelEllipsoid(radius_a, radius_b, radius_c);
  addJointAndBody(modelEllipsoid, jointModelEllipsoid, 0, pos, "ellipsoid", inertia);

  // Create Composite model (Tx, Ty, Tz, Rx, Ry, Rz)
  Model modelComposite;
  JointModelComposite jComposite;
  jComposite.addJoint(JointModelPX());
  jComposite.addJoint(JointModelPY());
  jComposite.addJoint(JointModelPZ());
  jComposite.addJoint(JointModelRX());
  jComposite.addJoint(JointModelRY());
  jComposite.addJoint(JointModelRZ());
  addJointAndBody(modelComposite, jComposite, 0, pos, "composite", inertia);

  Data dataEllipsoid(modelEllipsoid);
  Data dataComposite(modelComposite);

  // Test positions of ellispoid vs composite
  Eigen::VectorXd q_ellipsoid(3);
  q_ellipsoid << 1.0, 2.0, 3.0; // rx, ry, rz
  Eigen::Vector3d t = jointModelEllipsoid.computeTranslations(q_ellipsoid);
  Eigen::VectorXd q_composite(6);
  q_composite << t, q_ellipsoid;

  forwardKinematics(modelEllipsoid, dataEllipsoid, q_ellipsoid);
  forwardKinematics(modelComposite, dataComposite, q_composite);

  BOOST_CHECK(dataEllipsoid.oMi[1].isApprox(dataComposite.oMi[1]));

  // Velocity test
  Eigen::VectorXd qdot_ellipsoid(3);
  qdot_ellipsoid << 0.1, 0.2, 0.3;

  Eigen::Vector3d v_linear =
    jointModelEllipsoid.computeTranslationVelocities(q_ellipsoid, qdot_ellipsoid);
  Eigen::VectorXd qdot_composite(6);
  qdot_composite << v_linear, qdot_ellipsoid;

  forwardKinematics(modelEllipsoid, dataEllipsoid, q_ellipsoid, qdot_ellipsoid);
  forwardKinematics(modelComposite, dataComposite, q_composite, qdot_composite);

  BOOST_CHECK(dataEllipsoid.v[1].toVector().isApprox(dataComposite.v[1].toVector()));

  // Acceleration test
  Eigen::VectorXd qddot_ellipsoid(3);
  qddot_ellipsoid << 0.01, 0.02, 0.03;
  Eigen::Vector3d a_linear = jointModelEllipsoid.computeTranslationAccelerations(
    q_ellipsoid, qdot_ellipsoid, qddot_ellipsoid);
  Eigen::VectorXd qddot_composite(6);
  qddot_composite << a_linear, qddot_ellipsoid;

  forwardKinematics(modelEllipsoid, dataEllipsoid, q_ellipsoid, qdot_ellipsoid, qddot_ellipsoid);
  forwardKinematics(modelComposite, dataComposite, q_composite, qdot_composite, qddot_composite);

  BOOST_CHECK(dataEllipsoid.a[1].toVector().isApprox(dataComposite.a[1].toVector()));

  // Test RNEA - spatial forces and torques should match
  rnea(modelEllipsoid, dataEllipsoid, q_ellipsoid, qdot_ellipsoid, qddot_ellipsoid);
  rnea(modelComposite, dataComposite, q_composite, qdot_composite, qddot_composite);
  BOOST_CHECK(dataEllipsoid.f[1].isApprox(dataComposite.f[1]));

  // Need joint data to get both motion subspaces S_comp and S_ell
  JointDataComposite jdata_c = jComposite.createData();
  jComposite.setIndexes(0, 0, 0);
  jComposite.calc(jdata_c, q_composite, qdot_composite);

  JointDataEllipsoid jdata_e = jointModelEllipsoid.createData();
  jointModelEllipsoid.setIndexes(0, 0, 0);
  jointModelEllipsoid.calc(jdata_e, q_ellipsoid, qdot_ellipsoid);

  const Eigen::Matrix<double, 6, 6> S_comp = jdata_c.S.matrix(); // 6x6 (Tx,Ty,Tz,Rx,Ry,Rz)
  const Eigen::Matrix<double, 6, 3> S_ell = jdata_e.S.matrix();  // 6x3 (ellipsoid)

  // Compute the Jacobian mapping from ellipsoid generalized velocities to composite ones:
  // qdot_comp = J * qdot_ell with J = (S_comp^T * S_comp)^-1 * S_comp^T * S_ell
  Eigen::MatrixXd J = (S_comp.transpose() * S_comp).ldlt().solve(S_comp.transpose() * S_ell);

  // Now map the torques back (dual mapping)
  // τ_ell = J^T * τ_comp
  Eigen::VectorXd tau_proj = J.transpose() * dataComposite.tau;

  // Check the numerical match
  BOOST_CHECK_MESSAGE(
    dataEllipsoid.tau.isApprox(tau_proj, 1e-6),
    "Projected composite torques do not match ellipsoid torques.\n"
      << "Expected: " << dataEllipsoid.tau.transpose() << "\nGot: " << tau_proj.transpose());

  // Test ABA (Articulated-Body Algorithm)
  Eigen::VectorXd aAbaEllipsoid = aba(
    modelEllipsoid, dataEllipsoid, q_ellipsoid, qdot_ellipsoid, dataEllipsoid.tau,
    Convention::WORLD);
  BOOST_CHECK(aAbaEllipsoid.isApprox(qddot_ellipsoid));
}

BOOST_AUTO_TEST_SUITE_END()
