//
// Copyright (c) 2026 INRIA
// Copyright (c) 2026 ISIR
//

#pragma once

// IWYU pragma: private, include "pinocchio/algorithm/spline-regression.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/algorithm/spline-regression.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{

  ///
  /// \brief Structure containing all the settings parameters for the spline
  ///        regression algorithm (see fitSplineJoint).
  ///
  /// \tparam _Scalar Scalar type of the settings parameters.
  ///
  /// The Gauss-Newton iterations stop as soon as the data residual
  /// \f$\sum_i \|u_i\|^2\f$ falls below \ref absolute_accuracy, when its relative
  /// improvement between two iterates falls below \ref relative_accuracy, or
  /// after \ref max_iter iterations.
  ///
  /// Setting \ref smoothing_weight > 0 enables the Hessian-based smoother of
  /// Lee & Terzopoulos (2008), Sec. 6.3: the linear system blends the data-fit
  /// term with a penalty on the joint Hessian \f$\mathrm{d}S/\mathrm{d}q\f$,
  /// which damps spurious oscillations when fitting noisy data.
  ///
  template<typename _Scalar>
  struct SplineRegressionSettingsTpl
  {
    typedef _Scalar Scalar;

    /// \brief Default constructor.
    SplineRegressionSettingsTpl()
    : absolute_accuracy(Eigen::NumTraits<Scalar>::dummy_precision())
    , relative_accuracy(Eigen::NumTraits<Scalar>::dummy_precision())
    , damping(Scalar(0))
    , smoothing_weight(Scalar(0))
    , smoothing_drag(Scalar(1))
    , max_iter(100)
    , absolute_residual(Scalar(-1))
    , relative_residual(Scalar(-1))
    , iter(0)
    {
    }

    ///
    /// \brief Constructor with the main setting parameters.
    ///
    SplineRegressionSettingsTpl(const Scalar accuracy, const Scalar damping, const int max_iter)
    : absolute_accuracy(accuracy)
    , relative_accuracy(accuracy)
    , damping(damping)
    , smoothing_weight(Scalar(0))
    , smoothing_drag(Scalar(1))
    , max_iter(max_iter)
    , absolute_residual(Scalar(-1))
    , relative_residual(Scalar(-1))
    , iter(0)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(accuracy >= Scalar(0)) && "Accuracy must be positive.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        check_expression_if_real<Scalar>(damping >= Scalar(0)) && "Damping must be positive.");
      assert(max_iter >= 1 && "max_iter must be greater or equal to 1");
    }

    // parameters

    /// \brief Absolute accuracy: stop when the residual goes below this value.
    Scalar absolute_accuracy;

    /// \brief Relative accuracy: stop when the relative residual improvement
    ///        between two iterates goes below this value.
    Scalar relative_accuracy;

    /// \brief Tikhonov damping (\f$\lambda \geq 0\f$) added to the normal matrix.
    ///        Biases the control-frame corrections toward zero; useful on noisy
    ///        or sparse data.
    Scalar damping;

    /// \brief Hessian smoothing weight \f$w \in [0, 1)\f$ (0 disables smoothing).
    Scalar smoothing_weight;

    /// \brief Hessian smoothing drag rate \f$c \in (0, 1]\f$ (1 tries to cancel
    ///        the joint Hessian in a single Newton step).
    Scalar smoothing_drag;

    /// \brief Maximal number of Gauss-Newton iterations.
    int max_iter;

    // data that can be modified by the algorithm

    /// \brief Final data residual \f$\sum_i \|u_i\|^2\f$.
    Scalar absolute_residual;

    /// \brief Relative residual improvement at the last iterate.
    Scalar relative_residual;

    /// \brief Total number of iterations performed.
    int iter;
  };

  namespace internal
  {
    /// \brief Evaluate the cumulative B-spline placement
    ///        \f$X(q_i) = C_0 \prod_j \exp_6(v_j \tilde{B}_j(q_i))\f$ for the
    ///        row \p row of the cumulative basis matrix \p cumulative_basis.
    template<typename Scalar, int Options, typename MatrixLike>
    SE3Tpl<Scalar, Options> evaluateCumulativeSpline(
      const std::vector<SE3Tpl<Scalar, Options>> & ctrl_frames,
      const std::vector<MotionTpl<Scalar, Options>> & relative_motions,
      const Eigen::MatrixBase<MatrixLike> & cumulative_basis,
      const Eigen::Index row)
    {
      typedef SE3Tpl<Scalar, Options> SE3;
      SE3 placement = ctrl_frames[0];
      for (size_t j = 1; j < ctrl_frames.size(); ++j)
        placement =
          placement * exp6(relative_motions[j - 1] * cumulative_basis(row, Eigen::Index(j)));
      return placement;
    }
  } // namespace internal

  template<typename Scalar, int Options, typename ConfigVectorType>
  JointModelSplineTpl<Scalar, Options> fitSplineJoint(
    const std::vector<SE3Tpl<Scalar, Options>> & data_frames,
    const Eigen::MatrixBase<ConfigVectorType> & q_data,
    const size_t nbControlFrames,
    const size_t degree,
    SplineRegressionSettingsTpl<Scalar> & settings)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(ConfigVectorType);
    typedef SE3Tpl<Scalar, Options> SE3;
    typedef MotionTpl<Scalar, Options> Motion;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> MatrixXs;

    const Eigen::Index n_data = q_data.size();
    const Eigen::Index n_ctrl = Eigen::Index(nbControlFrames);

    PINOCCHIO_CHECK_ARGUMENT_SIZE(Eigen::Index(data_frames.size()), n_data);
    if (nbControlFrames <= degree)
      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "fitSplineJoint - Number of control frames must be greater than the spline degree.");
    if (n_data < n_ctrl)
      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "fitSplineJoint - Number of data frames must be greater or equal to the number of "
        "control frames.");
    for (Eigen::Index i = 1; i < n_data; ++i)
    {
      if (check_expression_if_real<Scalar>(q_data[i] <= q_data[i - 1]))
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "fitSplineJoint - q_data must be strictly increasing.");
    }
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      check_expression_if_real<Scalar>(
        settings.smoothing_weight >= Scalar(0) && settings.smoothing_weight < Scalar(1))
      && "smoothing_weight must lie in [0, 1).");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      check_expression_if_real<Scalar>(
        settings.smoothing_drag > Scalar(0) && settings.smoothing_drag <= Scalar(1))
      && "smoothing_drag must lie in (0, 1].");
    PINOCCHIO_CHECK_INPUT_ARGUMENT(
      check_expression_if_real<Scalar>(settings.damping >= Scalar(0))
      && "damping must be positive.");

    // Open uniform knots over the data parameter range: the fitted joint is
    // driven directly in the data parameter (e.g. a physical joint angle).
    const VectorXs knots =
      internal::generateOpenUniformKnots(q_data[0], q_data[n_data - 1], nbControlFrames, degree);

    const bool with_smoothing =
      check_expression_if_real<Scalar>(settings.smoothing_weight > Scalar(0));
    const Scalar weight = settings.smoothing_weight;
    const Scalar drag = settings.smoothing_drag;

    // Design matrices: basis values and, for the smoother, first and second
    // basis derivatives at the sample parameters. They only depend on the knot
    // vector, so they are computed once. The cumulative (right-summed) variants
    // feed the cumulative B-spline formula and its Hessian recursion.
    MatrixXs basis(MatrixXs::Zero(n_data, n_ctrl));
    MatrixXs basis_der2;
    MatrixXs cumulative_basis(n_data, n_ctrl);
    MatrixXs cumulative_der, cumulative_der2;
    if (with_smoothing)
    {
      basis_der2 = MatrixXs::Zero(n_data, n_ctrl);
      cumulative_der.resize(n_data, n_ctrl);
      cumulative_der2.resize(n_data, n_ctrl);
    }
    MatrixXs basis_der(with_smoothing ? n_data : 0, with_smoothing ? n_ctrl : 0);
    for (Eigen::Index i = 0; i < n_data; ++i)
    {
      const Scalar q_i = q_data[i];
      for (Eigen::Index j = 0; j < n_ctrl; ++j)
      {
        basis(i, j) = internal::bsplineBasis(size_t(j), degree, q_i, knots);
        if (with_smoothing)
        {
          basis_der(i, j) = internal::bsplineBasisDerivative(size_t(j), degree, q_i, knots);
          basis_der2(i, j) = internal::bsplineBasisDerivative2(size_t(j), degree, q_i, knots);
        }
      }
      // Cumulative sums from the right: tilde_B_j = sum_{l >= j} B_l.
      Scalar acc = Scalar(0), acc_der = Scalar(0), acc_der2 = Scalar(0);
      for (Eigen::Index j = n_ctrl - 1; j >= 0; --j)
      {
        acc += basis(i, j);
        cumulative_basis(i, j) = acc;
        if (with_smoothing)
        {
          acc_der += basis_der(i, j);
          cumulative_der(i, j) = acc_der;
          acc_der2 += basis_der2(i, j);
          cumulative_der2(i, j) = acc_der2;
        }
      }
    }

    // Normal matrix of the (simplified-Jacobian) Gauss-Newton step. It is
    // constant across iterations, so it is factorized once.
    MatrixXs normal_matrix(n_ctrl, n_ctrl);
    if (with_smoothing)
      normal_matrix.noalias() = (Scalar(1) - weight) * (basis.transpose() * basis)
                                + weight * (basis_der2.transpose() * basis_der2);
    else
      normal_matrix.noalias() = basis.transpose() * basis;
    normal_matrix.diagonal().array() += settings.damping;
    const Eigen::LDLT<MatrixXs> ldlt(normal_matrix);

    // Initialize each control frame with the data frame whose parameter is the
    // closest to the control frame Greville abscissa. This starts the iterations
    // near the solution and speeds up convergence w.r.t. an identity guess.
    std::vector<SE3> ctrl_frames(nbControlFrames);
    for (Eigen::Index k = 0; k < n_ctrl; ++k)
    {
      const Scalar greville =
        knots.segment(k + 1, Eigen::Index(degree)).sum() / Scalar(Eigen::Index(degree));
      Eigen::Index nearest = 0;
      (q_data.array() - greville).abs().minCoeff(&nearest);
      ctrl_frames[size_t(k)] = data_frames[size_t(nearest)];
    }

    // Gauss-Newton loop (Lee & Terzopoulos, Sec. 6.2-6.3, simplified Jacobian).
    std::vector<Motion> relative_motions(nbControlFrames - 1);
    MatrixXs residuals(n_data, 6);
    MatrixXs hessians(with_smoothing ? n_data : 0, with_smoothing ? 6 : 0);
    MatrixXs rhs(n_ctrl, 6), corrections(n_ctrl, 6);
    Scalar previous_cost = Scalar(-1);

    settings.relative_residual = Scalar(-1);
    for (int it = 0;; ++it)
    {
      for (size_t j = 1; j < nbControlFrames; ++j)
        relative_motions[j - 1] = log6(ctrl_frames[j - 1].inverse() * ctrl_frames[j]);

      for (Eigen::Index i = 0; i < n_data; ++i)
      {
        const SE3 placement =
          internal::evaluateCumulativeSpline(ctrl_frames, relative_motions, cumulative_basis, i);
        residuals.row(i) =
          log6(placement.inverse() * data_frames[size_t(i)]).toVector().transpose();
      }

      const Scalar cost = residuals.squaredNorm();
      settings.absolute_residual = cost;
      settings.iter = it;
      if (it > 0)
        settings.relative_residual =
          math::fabs(previous_cost - cost)
          / math::max(previous_cost, Eigen::NumTraits<Scalar>::dummy_precision());
      previous_cost = cost;

      if (check_expression_if_real<Scalar>(cost < settings.absolute_accuracy))
        break;
      if (
        it > 0
        && check_expression_if_real<Scalar>(
          settings.relative_residual < settings.relative_accuracy))
        break;
      if (it >= settings.max_iter)
        break;

      if (with_smoothing)
      {
        // Joint Hessian h(q_i) = dS/dq at the sample parameters, computed with
        // the same recursion as JointModelSplineTpl::computeTransformations.
        for (Eigen::Index i = 0; i < n_data; ++i)
        {
          Motion S(Motion::Zero()), S_der(Motion::Zero());
          for (size_t j = 1; j < nbControlFrames; ++j)
          {
            const Motion & v_j = relative_motions[j - 1];
            const SE3 factor(exp6(v_j * cumulative_basis(i, Eigen::Index(j))));
            const Motion S_der_next =
              v_j * cumulative_der2(i, Eigen::Index(j))
              + factor.actInv(S_der + S.cross(v_j) * cumulative_der(i, Eigen::Index(j)));
            S = v_j * cumulative_der(i, Eigen::Index(j)) + factor.actInv(S);
            S_der = S_der_next;
          }
          hessians.row(i) = S_der.toVector().transpose();
        }
        rhs.noalias() = (Scalar(1) - weight) * (basis.transpose() * residuals)
                        - (drag * weight) * (basis_der2.transpose() * hessians);
      }
      else
        rhs.noalias() = basis.transpose() * residuals;

      corrections = ldlt.solve(rhs);
      for (Eigen::Index k = 0; k < n_ctrl; ++k)
        ctrl_frames[size_t(k)] =
          ctrl_frames[size_t(k)] * exp6(Motion(corrections.row(k).transpose()));
    }

    return JointModelSplineTpl<Scalar, Options>(ctrl_frames, knots, degree);
  }

  template<typename Scalar, int Options, typename ConfigVectorType>
  JointModelSplineTpl<Scalar, Options> fitSplineJoint(
    const std::vector<SE3Tpl<Scalar, Options>> & data_frames,
    const Eigen::MatrixBase<ConfigVectorType> & q_data,
    const size_t nbControlFrames,
    const size_t degree)
  {
    SplineRegressionSettingsTpl<Scalar> settings;
    return fitSplineJoint(data_frames, q_data, nbControlFrames, degree, settings);
  }

} // namespace pinocchio
