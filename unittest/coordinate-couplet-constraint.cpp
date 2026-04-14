//
// Copyright (c) 2026 INRIA
//

#include "pinocchio/constraints.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include "constraints/jacobians-checker.hpp"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iterator>

using namespace pinocchio;

typedef CoordinateCoupletConstraintModel::EigenIndexVector EigenIndexVector;
typedef CoordinateCoupletConstraintModel::BooleanVector BooleanVector;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constraint_constructor_and_calc)
{
  Model model;
  buildModels::manipulator(model);

  Data data(model);
  Eigen::VectorXd q = neutral(model);
  q[0] = 0.3;
  q[1] = -0.2;
  data.q_in = q;

  CoordinateCoupletConstraintModel constraint_model(model, 0, 1);
  CoordinateCoupletConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  BOOST_CHECK_EQUAL(constraint_model.residualSize(), 1);
  BOOST_CHECK_CLOSE(constraint_data.constraint_residual[0], q[0] - q[1], 1e-12);
  BOOST_CHECK_EQUAL(constraint_model.getCoordinate1Id(), 0);
  BOOST_CHECK_EQUAL(constraint_model.getCoordinate2Id(), 1);
}

BOOST_AUTO_TEST_CASE(cast_and_variant)
{
  Model model;
  buildModels::manipulator(model);

  CoordinateCoupletConstraintModel cm(model, 0, 1);
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

  Data data(model);

  CoordinateCoupletConstraintModel constraint_model(model, 0, 1);
  CoordinateCoupletConstraintData constraint_data(constraint_model);

  const double eps_fd = 1e-8;
  for (int i = 0; i < 100; ++i)
  {
    const Eigen::VectorXd q0 = randomConfiguration(model);
    data.q_in = q0;
    constraint_model.calc(model, data, constraint_data);

    Eigen::MatrixXd jacobian_matrix(constraint_model.residualSize(), model.nv);
    constraint_model.jacobian(model, data, constraint_data, jacobian_matrix);

    Data data_fd(model);
    CoordinateCoupletConstraintData constraint_data_fd(constraint_model);
    Eigen::MatrixXd jacobian_matrix_fd(constraint_model.residualSize(), model.nv);

    for (Eigen::Index k = 0; k < model.nv; ++k)
    {
      Eigen::VectorXd v_eps = Eigen::VectorXd::Zero(model.nv);
      v_eps[k] = eps_fd;
      const Eigen::VectorXd q_plus = integrate(model, q0, v_eps);
      data_fd.q_in = q_plus;

      constraint_model.calc(model, data_fd, constraint_data_fd);
      jacobian_matrix_fd.col(k) =
        (constraint_data_fd.constraint_residual - constraint_data.constraint_residual) / eps_fd;
    }

    BOOST_CHECK(jacobian_matrix.isApprox(jacobian_matrix_fd, math::sqrt(eps_fd)));
  }

  const Eigen::VectorXd q = randomConfiguration(model);
  data.q_in = q;
  constraint_model.calc(model, data, constraint_data);
  check_jacobians_operations(model, data, constraint_model, constraint_data);
}

BOOST_AUTO_TEST_CASE(check_maps_and_row_support)
{
  Model model;
  buildModels::manipulator(model);

  Data data(model);
  const Eigen::VectorXd q = randomConfiguration(model);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  data.q_in = q;

  CoordinateCoupletConstraintModel constraint_model(model, 0, 1);
  CoordinateCoupletConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  Eigen::MatrixXd J(constraint_model.residualSize(), model.nv);
  constraint_model.jacobian(model, data, constraint_data, J);

  Eigen::VectorXd lambda = Eigen::VectorXd::Random(constraint_model.residualSize());
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
  constraint_model.mapConstraintForceToJointTorques(model, data, constraint_data, lambda, tau);
  BOOST_CHECK(tau.isApprox(J.transpose() * lambda));

  Eigen::VectorXd constraint_velocity(constraint_model.residualSize());
  constraint_model.mapJointMotionsToConstraintMotion(model, data, constraint_data, v, constraint_velocity);
  BOOST_CHECK(constraint_velocity.isApprox(J * v));

  BooleanVector row_sparsity_pattern;
  EigenIndexVector row_indexes;
  constraint_model.getRowSparsityPattern(model, data, constraint_data, 0, row_sparsity_pattern);
  constraint_model.getRowIndexes(model, data, constraint_data, 0, row_indexes);

  EigenIndexVector expected_indexes;
  const auto & indexes1 = model.span_indexes_vector[constraint_model.getJoint1Id()];
  const auto & indexes2 = model.span_indexes_vector[constraint_model.getJoint2Id()];
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
  data.q_in = q;

  CoordinateCoupletConstraintModel constraint_model(model, 0, 1);
  CoordinateCoupletConstraintData constraint_data(constraint_model);
  constraint_model.calc(model, data, constraint_data);

  std::vector<ConstraintModel> constraints;
  constraints.emplace_back(constraint_model);
  computeJointMinimalOrdering(model, data, constraints);

  const JointIndex joint1_id = constraint_model.getJoint1Id();
  const JointIndex joint2_id = constraint_model.getJoint2Id();
  BOOST_CHECK(data.constraints_supported_dim[joint1_id] > 0);
  BOOST_CHECK(data.constraints_supported_dim[joint2_id] > 0);
  BOOST_CHECK(data.joint_coupling_info(Eigen::Index(joint1_id), Eigen::Index(joint2_id)));

  BOOST_CHECK_THROW(
    constraint_model.appendCouplingConstraintInertias(
      model, data, constraint_data, Eigen::VectorXd::Ones(1), WorldFrameTag()),
    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
