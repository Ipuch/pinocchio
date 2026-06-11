//
// Copyright (c) 2026 INRIA
// Copyright (c) 2026 ISIR
//

#include "pinocchio/algorithm/spline-regression.hpp"

#include <boost/test/unit_test.hpp>

using namespace pinocchio;

/// @brief Evaluate the placement of a standalone spline joint at parameter q.
static SE3 evaluateJoint(const JointModelSpline & jmodel, const double q)
{
  JointModelSpline jm = jmodel;
  jm.setIndexes(0, 0, 0);
  JointDataSpline jdata = jm.createData();
  Eigen::VectorXd q_vec(1);
  q_vec[0] = q;
  jm.calc(jdata, q_vec);
  return jdata.M;
}

/// @brief Evaluate the joint Hessian dS/dq of a standalone spline joint at parameter q.
static Motion evaluateJointHessian(const JointModelSpline & jmodel, const double q)
{
  JointModelSpline jm = jmodel;
  jm.setIndexes(0, 0, 0);
  JointDataSpline jdata = jm.createData();
  Eigen::VectorXd q_vec(1), v_vec(1);
  q_vec[0] = q;
  v_vec[0] = 1.0; // with unit velocity, the bias c equals dS/dq.
  jm.calc(jdata, q_vec, v_vec);
  return jdata.c;
}

BOOST_AUTO_TEST_SUITE(SplineRegression)

/// The fitted joint must trace a parabola y = x^2 exactly (a cubic B-spline
/// represents a degree-2 polynomial without error). The joint is driven
/// directly in x as the knots span the data range.
BOOST_AUTO_TEST_CASE(fitParabola)
{
  const int n_data = 80;
  std::vector<SE3> data_frames;
  Eigen::VectorXd q_data(n_data);
  for (int i = 0; i < n_data; ++i)
  {
    const double x = -1.0 + 2.0 * double(i) / double(n_data - 1);
    q_data[i] = x;
    data_frames.push_back(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d(x, x * x, 0.)));
  }

  SplineRegressionSettings settings;
  const JointModelSpline jmodel =
    fitSplineJoint(data_frames, q_data, size_t(9), size_t(3), settings);

  BOOST_CHECK(settings.absolute_residual < 1e-12);
  BOOST_CHECK_EQUAL(jmodel.min_q, -1.0);
  BOOST_CHECK_EQUAL(jmodel.max_q, 1.0);

  for (int i = 0; i <= 200; ++i)
  {
    const double x = -1.0 + 2.0 * double(i) / 200.0;
    const SE3 placement = evaluateJoint(jmodel, x);
    BOOST_CHECK_SMALL(placement.translation().x() - x, 1e-10);
    BOOST_CHECK_SMALL(placement.translation().y() - x * x, 1e-10);
    BOOST_CHECK(placement.rotation().isIdentity(1e-10));
  }
}

/// Round trip: sample placements from a known spline joint and refit with the
/// same degree and number of control frames. The fitted joint must reproduce
/// the reference placements everywhere (also between the samples).
BOOST_AUTO_TEST_CASE(roundTrip)
{
  std::vector<SE3> ctrl_frames;
  for (int k = 0; k < 6; ++k)
  {
    const double s = double(k) / 5.0;
    const Eigen::Vector3d rotation_vector(0.2 * s, -0.1 * s, 1.5 * s);
    const Eigen::Vector3d translation(0.3 * s, 0.05 * std::sin(3.0 * s), -0.4 * s * s);
    ctrl_frames.push_back(SE3(exp3(rotation_vector), translation));
  }
  const JointModelSpline jmodel_ref = JointModelSplineBuilder()
                                        .withDegree(3)
                                        .withControlFrameVector(ctrl_frames)
                                        .withOpenUniformKnots(0.0, 1.0)
                                        .build();

  const int n_data = 40;
  std::vector<SE3> data_frames;
  Eigen::VectorXd q_data(n_data);
  for (int i = 0; i < n_data; ++i)
  {
    q_data[i] = double(i) / double(n_data - 1);
    data_frames.push_back(evaluateJoint(jmodel_ref, q_data[i]));
  }

  SplineRegressionSettings settings;
  const JointModelSpline jmodel_fit =
    fitSplineJoint(data_frames, q_data, ctrl_frames.size(), size_t(3), settings);

  BOOST_CHECK(settings.absolute_residual < 1e-12);
  for (int i = 0; i <= 200; ++i)
  {
    const double q = double(i) / 200.0;
    const SE3 placement_ref = evaluateJoint(jmodel_ref, q);
    const SE3 placement_fit = evaluateJoint(jmodel_fit, q);
    const Motion error = log6(placement_fit.inverse() * placement_ref);
    BOOST_CHECK_SMALL(error.toVector().norm(), 1e-6);
  }
}

/// A degree-1 spline with two control frames fitted on pure-rotation data is a
/// revolute joint: the fit must converge immediately and reproduce the data.
BOOST_AUTO_TEST_CASE(fitRevolute)
{
  const double angle_max = 1.2;
  const int n_data = 10;
  std::vector<SE3> data_frames;
  Eigen::VectorXd q_data(n_data);
  for (int i = 0; i < n_data; ++i)
  {
    const double angle = angle_max * double(i) / double(n_data - 1);
    q_data[i] = angle;
    data_frames.push_back(SE3(
      Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()).toRotationMatrix(),
      Eigen::Vector3d::Zero()));
  }

  SplineRegressionSettings settings;
  const JointModelSpline jmodel =
    fitSplineJoint(data_frames, q_data, size_t(2), size_t(1), settings);

  BOOST_CHECK(settings.absolute_residual < 1e-12);
  for (int i = 0; i <= 50; ++i)
  {
    const double angle = angle_max * double(i) / 50.0;
    const SE3 placement = evaluateJoint(jmodel, angle);
    const SE3 reference(
      Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()).toRotationMatrix(),
      Eigen::Vector3d::Zero());
    BOOST_CHECK_SMALL(log6(placement.inverse() * reference).toVector().norm(), 1e-10);
  }
}

/// The Hessian smoother (Lee & Terzopoulos Sec. 6.3) must reduce the joint
/// Hessian magnitude on noisy data while keeping a reasonable data fit.
BOOST_AUTO_TEST_CASE(hessianSmoothing)
{
  const int n_data = 40;
  std::vector<SE3> data_frames;
  Eigen::VectorXd q_data(n_data);
  unsigned int seed = 42;
  for (int i = 0; i < n_data; ++i)
  {
    const double x = double(i) / double(n_data - 1);
    // Deterministic pseudo-noise on a smooth translation curve.
    seed = seed * 1103515245u + 12345u;
    const double noise = 0.02 * (double(seed % 1000u) / 1000.0 - 0.5);
    q_data[i] = x;
    data_frames.push_back(
      SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d(x, 0.2 * std::sin(2.0 * x) + noise, 0.)));
  }

  SplineRegressionSettings settings_raw;
  const JointModelSpline jmodel_raw =
    fitSplineJoint(data_frames, q_data, size_t(12), size_t(3), settings_raw);

  SplineRegressionSettings settings_smooth;
  settings_smooth.smoothing_weight = 0.1;
  const JointModelSpline jmodel_smooth =
    fitSplineJoint(data_frames, q_data, size_t(12), size_t(3), settings_smooth);

  double hessian_raw = 0., hessian_smooth = 0.;
  for (int i = 0; i <= 100; ++i)
  {
    const double q = double(i) / 100.0;
    hessian_raw += evaluateJointHessian(jmodel_raw, q).toVector().squaredNorm();
    hessian_smooth += evaluateJointHessian(jmodel_smooth, q).toVector().squaredNorm();
  }
  BOOST_CHECK(hessian_smooth < hessian_raw);

  // The smoothed fit must remain close to the (noisy) data.
  for (int i = 0; i < n_data; ++i)
  {
    const SE3 placement = evaluateJoint(jmodel_smooth, q_data[i]);
    const Motion error = log6(placement.inverse() * data_frames[size_t(i)]);
    BOOST_CHECK_SMALL(error.toVector().norm(), 0.1);
  }
}

/// Settings outputs must be filled by the algorithm.
BOOST_AUTO_TEST_CASE(settingsOutputs)
{
  const int n_data = 20;
  std::vector<SE3> data_frames;
  Eigen::VectorXd q_data(n_data);
  for (int i = 0; i < n_data; ++i)
  {
    const double x = double(i) / double(n_data - 1);
    q_data[i] = x;
    data_frames.push_back(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d(x, x * x, 0.)));
  }

  SplineRegressionSettings settings(1e-10, 1e-8, 50);
  fitSplineJoint(data_frames, q_data, size_t(6), size_t(3), settings);
  BOOST_CHECK(settings.iter >= 0);
  BOOST_CHECK(settings.iter <= 50);
  BOOST_CHECK(settings.absolute_residual >= 0.);
}

/// Invalid inputs must throw.
BOOST_AUTO_TEST_CASE(invalidInputs)
{
  std::vector<SE3> data_frames(10, SE3::Identity());
  Eigen::VectorXd q_data = Eigen::VectorXd::LinSpaced(10, 0., 1.);

  // Not enough control frames w.r.t. the degree.
  BOOST_CHECK_THROW(
    fitSplineJoint(data_frames, q_data, size_t(3), size_t(3)), std::invalid_argument);
  // More control frames than data frames.
  BOOST_CHECK_THROW(
    fitSplineJoint(data_frames, q_data, size_t(11), size_t(3)), std::invalid_argument);
  // Non increasing q_data.
  Eigen::VectorXd q_bad = q_data;
  q_bad[5] = q_bad[4];
  BOOST_CHECK_THROW(
    fitSplineJoint(data_frames, q_bad, size_t(5), size_t(3)), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
