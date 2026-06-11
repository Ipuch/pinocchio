//
// Copyright (c) 2026 INRIA
// Copyright (c) 2026 ISIR
//
#pragma once

// IWYU pragma: begin_keep
#include <cstddef>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include "pinocchio/macros.hpp"
#include "pinocchio/context.hpp"
#include "pinocchio/spatial.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/utils/check.hpp"
// IWYU pragma: end_keep

namespace pinocchio
{
  template<typename Scalar>
  struct SplineRegressionSettingsTpl;
  typedef SplineRegressionSettingsTpl<context::Scalar> SplineRegressionSettings;

  ///
  /// \brief Fit the control frames of a spline joint to a set of data placements
  ///        (Gauss-Newton regression in \f$se(3)\f$).
  ///
  /// Given data placements \f$X_1, \dots, X_N\f$ sampled at parameters
  /// \f$q_1 < \dots < q_N\f$, this algorithm solves
  /// \f[
  ///   \underset{C_0, \dots, C_{m}}{\text{argmin}}
  ///   \sum_{i=1}^{N} \left\| \log_6\!\left(X(q_i)^{-1} X_i\right) \right\|^2
  /// \f]
  /// over the spline control frames \f$C_k\f$ by Gauss-Newton iterations in
  /// \f$se(3)\f$ with the simplified (Adjoint-free) Jacobian of Lee & Terzopoulos,
  /// "Spline Joints for Multibody Dynamics", SIGGRAPH 2008, Sec. 6.2.
  /// Optionally, the Hessian-based smoother of Sec. 6.3 penalizes
  /// \f$\|\mathrm{d}S/\mathrm{d}q\|^2\f$ to trade data fit against joint smoothness
  /// (useful on noisy data); see SplineRegressionSettingsTpl.
  ///
  /// The returned joint uses an open uniform knot vector spanning
  /// \f$[q_1, q_N]\f$, so it is directly driven by the data parameter (e.g. a
  /// physical joint angle).
  ///
  /// \tparam Scalar Scalar type.
  /// \tparam Options Eigen alignment options.
  /// \tparam ConfigVectorType Type of the sample parameter vector.
  ///
  /// \param[in] data_frames The data placements to fit (size N >= nbControlFrames).
  /// \param[in] q_data The sample parameters, strictly increasing (size N).
  /// \param[in] nbControlFrames Number of control frames m + 1 (> degree).
  /// \param[in] degree Degree of the B-spline basis functions.
  /// \param[in,out] settings The regression settings; on output it contains the
  ///                final residual and the number of iterations performed.
  ///
  /// \returns The fitted spline joint model.
  ///
  template<typename Scalar, int Options, typename ConfigVectorType>
  JointModelSplineTpl<Scalar, Options> fitSplineJoint(
    const std::vector<SE3Tpl<Scalar, Options>> & data_frames,
    const Eigen::MatrixBase<ConfigVectorType> & q_data,
    const size_t nbControlFrames,
    const size_t degree,
    SplineRegressionSettingsTpl<Scalar> & settings);

  ///
  /// \brief Fit the control frames of a spline joint to a set of data placements,
  ///        with default regression settings.
  ///
  /// \copydetails fitSplineJoint(const std::vector<SE3Tpl<Scalar, Options>> &, const
  /// Eigen::MatrixBase<ConfigVectorType> &, const size_t, const size_t,
  /// SplineRegressionSettingsTpl<Scalar> &)
  ///
  template<typename Scalar, int Options, typename ConfigVectorType>
  JointModelSplineTpl<Scalar, Options> fitSplineJoint(
    const std::vector<SE3Tpl<Scalar, Options>> & data_frames,
    const Eigen::MatrixBase<ConfigVectorType> & q_data,
    const size_t nbControlFrames,
    const size_t degree = 3);

} // namespace pinocchio

// IWYU pragma: begin_exports
#include "pinocchio/src/algorithm/spline-regression.hxx"
// IWYU pragma: end_exports
