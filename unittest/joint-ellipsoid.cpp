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

BOOST_AUTO_TEST_CASE(vsFreeFlyer)
{
  using namespace pinocchio;
  typedef SE3::Vector3 Vector3;
  typedef Eigen::Matrix<double, 6, 1> Vector6;
  typedef Eigen::Matrix<double, 7, 1> VectorFF;
  typedef SE3::Matrix3 Matrix3;

  Model modelEllipsoid, modelFreeflyer;

  Inertia inertia(1., Vector3(0.5, 0., 0.0), Matrix3::Identity());
  SE3 pos(1);
  pos.translation() = SE3::LinearType(1., 0., 0.);

  addJointAndBody(modelEllipsoid, JointModelEllipsoid(1, 2,3), 0, pos, "ellipsoid", inertia);
  addJointAndBody(modelFreeflyer, JointModelFreeFlyer(), 0, pos, "free-flyer", inertia);

  Data dataEllipsoid(modelEllipsoid);

  Eigen::VectorXd q = Eigen::VectorXd::Ones(modelEllipsoid.nq);
  Eigen::VectorXd v = Eigen::VectorXd::Ones(modelEllipsoid.nv);

  forwardKinematics(modelEllipsoid, dataEllipsoid, q, v);

  Eigen::VectorXd tauEllipsoid = Eigen::VectorXd::Ones(modelEllipsoid.nv);

  Eigen::VectorXd aEllipsoid = Eigen::VectorXd::Ones(modelEllipsoid.nv);

  // ForwardDynamics == aba
  Eigen::VectorXd aAbaEllipsoid =
  aba(modelEllipsoid, dataEllipsoid, q, v, tauEllipsoid, Convention::WORLD);

  // Calculer jdata.S().transpose() * data.f[i]
}

// BOOST_AUTO_TEST_CASE(vsSphericalZYX)
// {
//   using namespace pinocchio;
//   typedef SE3::Vector3 Vector3;
//   typedef Eigen::Matrix<double, 6, 1> Vector6;
//   typedef SE3::Matrix3 Matrix3;

//   // Build models using ModelGraph to enable configuration conversion
//   graph::ModelGraph g;
//   Inertia inertia(1., Vector3(0.5, 0., 0.0), Matrix3::Identity());
//   SE3 pos(1);
//   pos.translation() = SE3::LinearType(1., 0., 0.);

//   g.addBody("root_body", inertia);
//   g.addBody("end_body", inertia);

//   // Create SphericalZYX joint (uses ZYX Euler angles)
//   g.addJoint(
//     "spherical_joint",
//     graph::JointSphericalZYX(),
//     "root_body",
//     SE3::Identity(),
//     "end_body",
//     pos);

//   // Build SphericalZYX model
//   const auto forward_build = graph::buildModelWithBuildInfo(g, "root_body", SE3::Identity());
//   const Model & modelSphericalZYX = forward_build.model;
//   Data dataSphericalZYX(modelSphericalZYX);

//   Eigen::VectorXd q_zyx(3);
//   q_zyx << 0.5, 1.2, -0.8;
//   Eigen::VectorXd v_zyx(3);
//   v_zyx << 0.1, -0.3, 0.7;
//   forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx);
//   updateFramePlacements(modelSphericalZYX, dataSphericalZYX);
  
//   auto end_index = modelSphericalZYX.getFrameId("end_body", BODY);
//   auto X_end = dataSphericalZYX.oMf[end_index];
//   const auto backward_build =
//     graph::buildModelWithBuildInfo(g, "end_body", X_end);

//   const Model & XYZlacirehpSledom = backward_build.model;

//   // Create converter from SphericalZYX to Ellipsoid
//   auto converter = graph::createConverter(
//     modelSphericalZYX, XYZlacirehpSledom, forward_build.build_info, backward_build.build_info);

//   // Convert to Ellipsoid configuration
//   Eigen::VectorXd q_xyz = Eigen::VectorXd::Zero(XYZlacirehpSledom.nq);
//   Eigen::VectorXd v_xyz = Eigen::VectorXd::Zero(XYZlacirehpSledom.nv);
//   converter.convertConfigurationVector(q_zyx, q_xyz);
//   converter.convertTangentVector(q_zyx, v_zyx, v_xyz);

//   std::cout << "\n=== Configuration Conversion ===" << std::endl;
//   std::cout << "q_zyx (SphericalZYX): " << q_zyx.transpose() << std::endl;
//   std::cout << "q_xyz (converted): " << q_xyz.transpose() << std::endl;
//   std::cout << "v_zyx (SphericalZYX): " << v_zyx.transpose() << std::endl;
//   std::cout << "v_xyz (converted): " << v_xyz.transpose() << std::endl;
  
//   Model modelEllipsoid;
//   addJointAndBody(modelEllipsoid, JointModelEllipsoid(0, 0, 0), 0, pos, "ellipsoid", inertia);
//   Data dataEllipsoid(modelEllipsoid);

//   // Test forward kinematics with converted configurations
//   forwardKinematics(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz);
//   forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx);

//   // Check that the transformations are identical
//   std::cout << "\n=== Forward Kinematics Comparison ===" << std::endl;
//   std::cout << "oMi (Ellipsoid):\n" << dataEllipsoid.oMi[1] << std::endl;
//   std::cout << "oMi (SphericalZYX):\n" << dataSphericalZYX.oMi[1] << std::endl;
//   BOOST_CHECK(dataEllipsoid.oMi[1].isApprox(dataSphericalZYX.oMi[1]));
  
//   BOOST_CHECK(dataEllipsoid.liMi[1].isApprox(dataSphericalZYX.liMi[1]));
  
//   // Check that velocities match
//   BOOST_CHECK(dataEllipsoid.v[1].toVector().isApprox(dataSphericalZYX.v[1].toVector()));

//   // Test computeAllTerms
//   computeAllTerms(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz);
//   computeAllTerms(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx);

//   BOOST_CHECK(dataEllipsoid.Ycrb[1].matrix().isApprox(dataSphericalZYX.Ycrb[1].matrix()));
  
//   BOOST_CHECK(dataEllipsoid.f[1].toVector().isApprox(dataSphericalZYX.f[1].toVector()));
  
//   BOOST_CHECK(dataEllipsoid.nle.isApprox(dataSphericalZYX.nle));
  
//   BOOST_CHECK(dataEllipsoid.com[0].isApprox(dataSphericalZYX.com[0]));

//   // Test inverse dynamics (RNEA)
//   Eigen::VectorXd aEllipsoid = Eigen::VectorXd::Ones(modelEllipsoid.nv);
//   Eigen::VectorXd aSphericalZYX = Eigen::VectorXd::Ones(modelSphericalZYX.nv);
  
//   Eigen::VectorXd tauEllipsoid = rnea(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz, aEllipsoid);
//   Eigen::VectorXd tauSphericalZYX = rnea(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx, aSphericalZYX);
  
//   BOOST_CHECK(tauEllipsoid.isApprox(tauSphericalZYX));

//   // Test forward dynamics (ABA)
//   Eigen::VectorXd tau = Eigen::VectorXd::Ones(modelEllipsoid.nv);
  
//   Eigen::VectorXd aAbaEllipsoid = aba(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz, tau, Convention::WORLD);
//   Eigen::VectorXd aAbaSphericalZYX = aba(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx, tau, Convention::WORLD);
  
//   BOOST_CHECK(aAbaEllipsoid.isApprox(aAbaSphericalZYX));

//   // Test with LOCAL convention
//   aAbaEllipsoid = aba(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz, tau, Convention::LOCAL);
//   aAbaSphericalZYX = aba(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx, tau, Convention::LOCAL);
  
//   BOOST_CHECK(aAbaEllipsoid.isApprox(aAbaSphericalZYX));

//   // Test with different configurations
//   q_zyx << 0.2, -0.5, 1.1;
//   v_zyx << 0.3, 0.1, -0.2;
  
//   converter.convertConfigurationVector(q_zyx, q_xyz);
//   converter.convertTangentVector(q_zyx, v_zyx, v_xyz);
  
//   forwardKinematics(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz);
//   forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx);

//   BOOST_CHECK(dataEllipsoid.oMi[1].isApprox(dataSphericalZYX.oMi[1]));
  
//   BOOST_CHECK(dataEllipsoid.v[1].toVector().isApprox(dataSphericalZYX.v[1].toVector()));

//   tauEllipsoid = rnea(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz, aEllipsoid);
//   tauSphericalZYX = rnea(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx, aSphericalZYX);
  
//   BOOST_CHECK(tauEllipsoid.isApprox(tauSphericalZYX));

//   aAbaEllipsoid = aba(modelEllipsoid, dataEllipsoid, q_xyz, v_xyz, tau, Convention::WORLD);
//   aAbaSphericalZYX = aba(modelSphericalZYX, dataSphericalZYX, q_zyx, v_zyx, tau, Convention::WORLD);
  
//   BOOST_CHECK(aAbaEllipsoid.isApprox(aAbaSphericalZYX));
  
// }

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
  q_s << 0.5, 1.2, -0.8;  // Z=0.5, Y=1.2, X=-0.8
  Eigen::VectorXd qd_s(3);
  qd_s << 0.1, -0.3, 0.7;  // ZYX angle velocities
  
  // Compute the rotation matrix from ZYX angles
  forwardKinematics(modelSphericalZYX, dataSphericalZYX, q_s, qd_s);
  const Matrix3 & R = dataSphericalZYX.oMi[1].rotation();
  const Motion & spatial_vel_zyx = dataSphericalZYX.v[1];

  std::cout << "\n=== Target Rotation Matrix (from ZYX) ===" << std::endl;
  std::cout << R << std::endl;
  
  // Extract XYZ Euler angles from the rotation matrix
  // For XYZ convention: R = Rx(x) * Ry(y) * Rz(z)
  // We need to solve for x, y, z
  
  Eigen::Vector3d q_e;
  
  // From the Ellipsoid rotation matrix structure:
  // R(0,2) = sin(y)
  // R(1,2) = -sin(x)*cos(y)
  // R(2,2) = cos(x)*cos(y)
  
  double sy = R(0, 2);
  q_e(1) = std::asin(sy);  // y angle
  
  double cy = std::cos(q_e(1));
  
  if (std::abs(cy) > 1e-6) {
    // Not at singularity
    q_e(0) = std::atan2(-R(1, 2) / cy, R(2, 2) / cy);  // x angle
    q_e(2) = std::atan2(-R(0, 1) / cy, R(0, 0) / cy);  // z angle
  } else {
    // Gimbal lock - choose x = 0
    q_e(0) = 0.0;
    q_e(2) = std::atan2(R(1, 0), R(1, 1));
  }
  
  std::cout << "\n=== Configuration ===" << std::endl;
  std::cout << "q_s (Z,Y,X): " << q_s.transpose() << std::endl;
  std::cout << "q_e (X,Y,Z): " << q_e.transpose() << std::endl;

  // Get the motion subspace matrices (which give us the Jacobians)
  JointModelSphericalZYX jmodel_s;
  JointDataSphericalZYX jdata_s = jmodel_s.createData();
  jmodel_s.calc(jdata_s, q_s);
  
  JointModelEllipsoid jmodel_e(0, 0, 0);
  JointDataEllipsoid jdata_e = jmodel_e.createData();
  jmodel_e.calc(jdata_e, q_e);


  // The motion subspace S gives us: omega = S * v
  Matrix3 S_s = jdata_s.S.matrix().bottomRows<3>();  // Angular part
  Matrix3 S_e = jdata_e.S.matrix().bottomRows<3>();  // Angular part
  
  Eigen::Vector3d qd_e = S_e.inverse() * S_s * qd_s;

  Eigen::Vector3d w_s = S_s * qd_s;
  Eigen::Vector3d w_e = S_e * qd_e;
  

  std::cout << "\n=== Joint Velocities ===" << std::endl;
  std::cout << "qd_s (joint velocities): " << qd_s.transpose() << std::endl;
  std::cout << "qd_e (joint velocities): " << qd_e.transpose() << std::endl;
  
  std::cout << "\n=== Motion Subspace Matrices ===" << std::endl;
  std::cout << "S_zyx (full 6x3):\n" << jdata_s.S.matrix() << std::endl;
  std::cout << "\nS_xyz (full 6x3):\n" << jdata_e.S.matrix() << std::endl;
  
  std::cout << "\n=== Angular Velocities from Joint Velocities ===" << std::endl;
  std::cout << "omega from ZYX (S_s * qd_s): " << w_s.transpose() << std::endl;
  std::cout << "omega from XYZ (S_e * qd_e): " << w_e.transpose() << std::endl;
  
  BOOST_CHECK(w_s.isApprox(w_e));
  std::cout << "✓ Angular velocities from joint velocities match" << std::endl;

  // Compute forward kinematics with the converted configurations
  forwardKinematics(modelEllipsoid, dataEllipsoid, q_e, qd_e);
  
  // Also get the joint data to see body-frame velocities
  JointDataEllipsoid jdata_e_fk = jmodel_e.createData();
  JointDataEllipsoid jdata_e_fk2 = jmodel_e.createData();
  jmodel_e.calc(jdata_e_fk, q_e, qd_e);
  jmodel_e.calc(jdata_e_fk2, q_e);

  std::cout << "\n=== Motion Subspace Matrices after FK ===" << std::endl;
  std::cout << "S_e calc q:\n" << jdata_e_fk.S.matrix() << std::endl;
  std::cout << "\nS_e calc q v:\n" << jdata_e_fk2.S.matrix() << std::endl;
  
  BOOST_CHECK(jdata_e_fk.S.matrix().isApprox(jdata_e_fk2.S.matrix()));
  std::cout << "✓ Motion subspace matrices match" << std::endl;
  
  JointDataSphericalZYX jdata_s_fk = jmodel_s.createData();
  jmodel_s.calc(jdata_s_fk, q_s, qd_s);
  
  std::cout << "\n=== Joint-frame velocities (S * v) ===" << std::endl;
  Eigen::Matrix<double, 6, 1> joint_vel_e = jdata_e_fk.S.matrix() * qd_e;
  std::cout << "joint_vel_e: " << joint_vel_e.transpose() << std::endl;
  std::cout << "jdata_e_fk.v : " << jdata_e_fk.v.toVector().transpose() << std::endl;

  

  Eigen::Matrix<double, 6, 1> joint_vel_s = jdata_s_fk.S.matrix() * qd_s;
  std::cout << "Ellipsoid (S_xyz * qd_e): " << joint_vel_e.transpose() << std::endl;
  std::cout << "SphericalZYX (S_zyx * qd_s): " << joint_vel_s.transpose() << std::endl;

  std::cout << "\n=== Rotation Matrices ===" << std::endl;
  std::cout << "Ellipsoid rotation:\n" << dataEllipsoid.oMi[1].rotation() << std::endl;
  std::cout << "\nSphericalZYX rotation:\n" << dataSphericalZYX.oMi[1].rotation() << std::endl;
  
  std::cout << "\n=== Translations ===" << std::endl;
  std::cout << "Ellipsoid translation: " << dataEllipsoid.oMi[1].translation().transpose() << std::endl;
  std::cout << "SphericalZYX translation: " << dataSphericalZYX.oMi[1].translation().transpose() << std::endl;
  
  std::cout << "\n=== Spatial Velocities ===" << std::endl;
  std::cout << "Ellipsoid v:\n" << dataEllipsoid.v[1].toVector().transpose() << std::endl;
  std::cout << "SphericalZYX v:\n" << dataSphericalZYX.v[1].toVector().transpose() << std::endl;

  BOOST_CHECK(dataEllipsoid.v[1].toVector().isApprox(dataSphericalZYX.v[1].toVector()));
  std::cout << "✓ Spatial velocities match" << std::endl;
  
  BOOST_CHECK(dataEllipsoid.oMi[1].isApprox(dataSphericalZYX.oMi[1]));
  std::cout << "✓ Full oMi[1] matches" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
