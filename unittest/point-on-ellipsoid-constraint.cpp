//
// Copyright (c) 2026 INRIA
//

#include "pinocchio/constraints.hpp"

#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include "constraints/jacobians-checker.hpp"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iterator>

using namespace pinocchio;

typedef PointOnEllipsoidConstraintModel::EigenIndexVector EigenIndexVector;
typedef PointOnEllipsoidConstraintModel::BooleanVector BooleanVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_constructor_and_calc)
{
  Model model;
  buildModels::manipulator(model);

  Data data(model);
  const Eigen::VectorXd q = randomConfiguration(model);
  forwardKinematics(model, data, q);

  const JointIndex joint1_id = model.getJointId("shoulder2_joint");
  const JointIndex joint2_id = model.getJointId("wrist2_joint");
  const SE3 joint1_placement = SE3::Random();
  const SE3 joint2_placement = SE3::Random();
  Eigen::Vector3d semi_axes;
  semi_axes << 0.3, 0.4, 0.5;

  PointOnEllipsoidConstraintModel constraint_model(
    model, joint1_id, joint1_placement, joint2_id, joint2_placement, semi_axes);
  PointOnEllipsoidConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  BOOST_CHECK_EQUAL(constraint_model.residualSize(), 1);
  BOOST_CHECK(constraint_model.getSemiAxes().isApprox(semi_axes));
  BOOST_CHECK(constraint_model.joint1_id == joint1_id);
  BOOST_CHECK(constraint_model.joint2_id == joint2_id);

  const SE3 oMc1 = data.oMi[joint1_id] * joint1_placement;
  const SE3 oMc2 = data.oMi[joint2_id] * joint2_placement;
  const Eigen::Vector3d d_local = oMc2.rotation().transpose() * (oMc1.translation() - oMc2.translation());
  Eigen::Vector3d inv_sq;
  inv_sq << 1.0 / (semi_axes[0] * semi_axes[0]), 1.0 / (semi_axes[1] * semi_axes[1]),
    1.0 / (semi_axes[2] * semi_axes[2]);
  const double residual_ref = d_local.dot(inv_sq.cwiseProduct(d_local)) - 1.0;

  BOOST_CHECK_CLOSE(constraint_data.constraint_residual[0], residual_ref, 1e-10);
}

BOOST_AUTO_TEST_CASE(cast_and_variant)
{
  Model model;
  buildModels::manipulator(model);

  Eigen::Vector3d semi_axes;
  semi_axes << 0.3, 0.4, 0.5;
  PointOnEllipsoidConstraintModel cm(
    model, model.getJointId("shoulder2_joint"), SE3::Random(),
    model.getJointId("wrist2_joint"), SE3::Random(), semi_axes);

  const auto cm_cast_double = cm.cast<double>();
  BOOST_CHECK(cm_cast_double == cm);

  ConstraintModel generic_model(cm);
  ConstraintData generic_data = generic_model.createData();
  BOOST_CHECK_EQUAL(generic_model.residualSize(), 1);
  BOOST_CHECK(ConstraintData(cm.createData()) == generic_data);
}

BOOST_AUTO_TEST_CASE(constraint_jacobian)
{
  Model model;
  buildModels::manipulator(model);

  Eigen::Vector3d semi_axes;
  semi_axes << 0.3, 0.4, 0.5;

  PointOnEllipsoidConstraintModel constraint_model(
    model, model.getJointId("shoulder2_joint"), SE3::Random(),
    model.getJointId("wrist2_joint"), SE3::Random(), semi_axes);
  PointOnEllipsoidConstraintData constraint_data(constraint_model);

  const double eps_fd = 1e-8;
  for (int i = 0; i < 20; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model);

    Data data(model);
    forwardKinematics(model, data, q0);
    computeJointJacobians(model, data, q0);
    constraint_model.calc(model, data, constraint_data);

    Eigen::MatrixXd jacobian_matrix(constraint_model.residualSize(), model.nv);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

    Eigen::MatrixXd jacobian_matrix_fd(constraint_model.residualSize(), model.nv);

    for (Eigen::Index k = 0; k < model.nv; ++k)
    {
      Eigen::VectorXd v_eps = Eigen::VectorXd::Zero(model.nv);
      v_eps[k] = eps_fd;
      const Eigen::VectorXd q_plus = integrate(model, q0, v_eps);

      Data data_fd(model);
      PointOnEllipsoidConstraintData constraint_data_fd(constraint_model);
      forwardKinematics(model, data_fd, q_plus);
      constraint_model.calc(model, data_fd, constraint_data_fd);

      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
    }

    BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }

  const Eigen::VectorXd q = randomConfiguration(model);
  Data data(model);
  forwardKinematics(model, data, q);
  computeJointJacobians(model, data, q);
  constraint_model.calc(model, data, constraint_data);
  check_jacobians_operations(model, data, constraint_model, constraint_data);
}

BOOST_AUTO_TEST_CASE(check_maps_and_row_support)
{
  Model model;
  buildModels::manipulator(model);

  Data data(model);
  const Eigen::VectorXd q = randomConfiguration(model);
  forwardKinematics(model, data, q);
  computeJointJacobians(model, data, q);

  Eigen::Vector3d semi_axes;
  semi_axes << 0.3, 0.4, 0.5;
  PointOnEllipsoidConstraintModel constraint_model(
    model, model.getJointId("shoulder2_joint"), SE3::Random(),
    model.getJointId("wrist2_joint"), SE3::Random(), semi_axes);
  PointOnEllipsoidConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  Eigen::MatrixXd J(constraint_model.residualSize(), model.nv);
  constraint_model.jacobian(model, data, constraint_data, J);

  BooleanVector row_sparsity_pattern;
  EigenIndexVector row_indexes;
  constraint_model.getRowSparsityPattern(model, data, constraint_data, 0, row_sparsity_pattern);
  constraint_model.getRowIndexes(model, data, constraint_data, 0, row_indexes);

  EigenIndexVector expected_indexes;
  const auto & indexes1 = model.span_indexes_vector[constraint_model.joint1_id];
  const auto & indexes2 = model.span_indexes_vector[constraint_model.joint2_id];
  std::set_union(
    indexes1.begin(), indexes1.end(), indexes2.begin(), indexes2.end(),
    std::back_inserter(expected_indexes));

  BOOST_CHECK(row_indexes == expected_indexes);
  for (const Eigen::Index index : expected_indexes)
  {
    BOOST_CHECK(row_sparsity_pattern[index]);
  }
}

BOOST_AUTO_TEST_CASE(ordering_and_inertia_contract)
{
  Model model;
  buildModels::manipulator(model);

  Data data(model);
  const Eigen::VectorXd q = randomConfiguration(model);
  forwardKinematics(model, data, q);
  computeJointJacobians(model, data, q);

  Eigen::Vector3d semi_axes;
  semi_axes << 0.3, 0.4, 0.5;
  PointOnEllipsoidConstraintModel constraint_model(
    model, model.getJointId("shoulder2_joint"), SE3::Random(),
    model.getJointId("wrist2_joint"), SE3::Random(), semi_axes);
  PointOnEllipsoidConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  std::vector<ConstraintModel> constraints;
  constraints.emplace_back(constraint_model);
  computeJointMinimalOrdering(model, data, constraints);

  const JointIndex joint1_id = constraint_model.joint1_id;
  const JointIndex joint2_id = constraint_model.joint2_id;
  BOOST_CHECK(data.constraints_supported_dim[joint1_id] > 0);
  BOOST_CHECK(data.constraints_supported_dim[joint2_id] > 0);
  BOOST_CHECK(data.joint_coupling_info(Eigen::Index(joint1_id), Eigen::Index(joint2_id)));

  BOOST_CHECK_THROW(
    constraint_model.appendCouplingConstraintInertias(
      model, data, constraint_data, Eigen::VectorXd::Ones(1), WorldFrameTag()),
    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
