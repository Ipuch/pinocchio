//
// Copyright (c) 2026 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/constraints.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/constraints.hpp"
#endif // PINOCCHIO_LSP

#include <stdexcept>

namespace pinocchio
{

  // --------------------------------------------------------------
  // Cast
  // --------------------------------------------------------------
  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, PointOnEllipsoidConstraintModelTpl<Scalar, Options>>
  {
    typedef PointOnEllipsoidConstraintModelTpl<NewScalar, Options> type;
  };

  // --------------------------------------------------------------
  // Traits
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct traits<PointOnEllipsoidConstraintModelTpl<_Scalar, _Options>>
  {
    typedef PointOnEllipsoidConstraintModelTpl<_Scalar, _Options> ConstraintModel;
    typedef PointOnEllipsoidConstraintDataTpl<_Scalar, _Options> ConstraintData;

    typedef ConstraintModel Model;
    typedef ConstraintData Data;

    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    static constexpr ConstraintSizeType constraint_size_type = ConstraintSizeType::STATIC;

    static constexpr bool has_baumgarte_corrector = true;
    static constexpr bool has_set = true;
    static constexpr bool is_inequality_constraint = false;

    typedef FullSpaceConeTpl<Scalar, Options> ConstraintSet;
    typedef ZeroConeJordanOperationTpl<Scalar, Options> JordanOperation;
    typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;

    static constexpr int Size = 1;
    static constexpr int SymmetricConeSize = JordanOperation::ConeSize;
    static constexpr int SymmetricConeScalingSize = JordanOperation::ConeScalingSize;

    typedef Eigen::Matrix<Scalar, Size, 1, Options> ResidualVectorType;
    typedef Eigen::Matrix<Scalar, Size, Eigen::Dynamic, Eigen::RowMajor> JacobianMatrixType;
    typedef Eigen::Matrix<Scalar, SymmetricConeSize, 1, Options> ConeVectorType;
    typedef Eigen::Matrix<Scalar, SymmetricConeScalingSize, 1, Options> ConeScalingVectorType;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> VectorXs;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXs;

    template<typename InputMatrix>
    struct JacobianMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Size,
        InputMatrixPlain::ColsAtCompileTime,
        Eigen::RowMajor>
        type;
    };

    template<typename InputMatrix>
    struct JacobianTransposeMatrixProductReturnType
    {
      typedef typename InputMatrix::Scalar Scalar;
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(InputMatrix) InputMatrixPlain;
      typedef Eigen::Matrix<
        Scalar,
        Eigen::Dynamic,
        InputMatrixPlain::ColsAtCompileTime,
        InputMatrixPlain::Options>
        type;
    };
  };

  template<typename _Scalar, int _Options>
  struct traits<PointOnEllipsoidConstraintDataTpl<_Scalar, _Options>>
  : traits<PointOnEllipsoidConstraintModelTpl<_Scalar, _Options>>
  {
  };

  ///
  /// \brief Constraint enforcing that a point attached to body 1 lies on an
  ///        ellipsoid attached to body 2.
  ///
  /// The constraint is the scalar function \f$ c(q) = d^\top D d - 1 \f$ where
  /// \f$ d = {}^{c_2} R_o (p_1 - p_2) \f$ is the position of the point expressed
  /// in the ellipsoid local frame, and \f$ D = \mathrm{diag}(1/a^2, 1/b^2, 1/c^2) \f$
  /// is built from the three semi-axes \f$ (a, b, c) \f$ of the ellipsoid.
  ///
  /// \pre Before calling calc/jacobian, the user must have updated `data.oMi`
  /// (e.g. via `forwardKinematics`) and `data.J` (via `computeJointJacobians`).
  ///
  template<typename _Scalar, int _Options>
  struct PointOnEllipsoidConstraintModelTpl
  : BinaryKinematicsConstraintModelBase<PointOnEllipsoidConstraintModelTpl<_Scalar, _Options>>
  {
    typedef PointOnEllipsoidConstraintModelTpl Self;
    typedef BinaryKinematicsConstraintModelBase<Self> Base;
    typedef typename Base::BaseCommonParameters BaseCommonParameters;
    typedef ConstraintModelBase<Self> RootBase;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    static constexpr ConstraintSizeType constraint_size_type = traits<Self>::constraint_size_type;
    static constexpr bool has_baumgarte_corrector = traits<Self>::has_baumgarte_corrector;

    typedef typename traits<Self>::ConstraintSet ConstraintSet;
    typedef typename traits<Self>::JordanOperation JordanOperation;
    typedef typename traits<Self>::BaumgarteCorrectorParameters BaumgarteCorrectorParameters;

    static constexpr int Size = traits<Self>::Size;

    typedef typename traits<Self>::ResidualVectorType ResidualVectorType;
    typedef typename traits<Self>::JacobianMatrixType JacobianMatrixType;

    typedef typename Base::SE3 SE3;
    typedef typename Base::Vector3 Vector3;
    typedef typename Base::Matrix3 Matrix3;

    template<typename NewScalar, int NewOptions>
    friend struct PointOnEllipsoidConstraintModelTpl;

    template<typename NewScalar, int NewOptions>
    friend struct PointOnEllipsoidConstraintDataTpl;

    using RootBase::classname;
    using RootBase::residualSize;
    using typename RootBase::BooleanVector;
    using typename RootBase::EigenIndexVector;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    BaseCommonParameters & base_common_parameters()
    {
      return static_cast<BaseCommonParameters &>(*this);
    }

    const BaseCommonParameters & base_common_parameters() const
    {
      return static_cast<const BaseCommonParameters &>(*this);
    }

    /// \brief Default constructor. Semi-axes are set to one.
    PointOnEllipsoidConstraintModelTpl()
    : Base()
    , m_semi_axes(Vector3::Ones())
    , m_inv_squared_axes(Vector3::Ones())
    {
    }

    /// \brief Full constructor.
    /// \param[in] model Kinematic tree.
    /// \param[in] joint1_id Joint id supporting the point.
    /// \param[in] joint1_placement Placement of the point w.r.t. joint1.
    /// \param[in] joint2_id Joint id supporting the ellipsoid.
    /// \param[in] joint2_placement Placement of the ellipsoid frame w.r.t. joint2.
    /// \param[in] semi_axes Strictly positive semi-axes (a, b, c) of the ellipsoid,
    ///            expressed in the ellipsoid frame.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointOnEllipsoidConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const SE3 & joint1_placement,
      const JointIndex joint2_id,
      const SE3 & joint2_placement,
      const Vector3 & semi_axes)
    : Base(model, joint1_id, joint1_placement, joint2_id, joint2_placement)
    , m_semi_axes(semi_axes)
    {
      check_and_init_axes();
    }

    /// \brief Constructor with identity placements.
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    PointOnEllipsoidConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const JointIndex joint1_id,
      const JointIndex joint2_id,
      const Vector3 & semi_axes)
    : Base(model, joint1_id, joint2_id)
    , m_semi_axes(semi_axes)
    {
      check_and_init_axes();
    }

    template<typename NewScalar>
    typename CastType<NewScalar, PointOnEllipsoidConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, PointOnEllipsoidConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::template cast<NewScalar>(res);
      res.m_semi_axes = m_semi_axes.template cast<NewScalar>();
      res.m_inv_squared_axes = m_inv_squared_axes.template cast<NewScalar>();
      return res;
    }

    bool operator==(const PointOnEllipsoidConstraintModelTpl & other) const
    {
      return base() == other.base() && m_semi_axes == other.m_semi_axes;
    }

    bool operator!=(const PointOnEllipsoidConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    static std::string classnameImpl()
    {
      return std::string("PointOnEllipsoidConstraintModel");
    }

    std::string shortnameImpl() const
    {
      return classname();
    }

    ConstraintData createDataImpl() const
    {
      return ConstraintData(*this);
    }

    ConstraintSet setImpl(const ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      return ConstraintSet();
    }

    /// \brief Semi-axes (a, b, c) of the ellipsoid in its local frame.
    const Vector3 & getSemiAxes() const
    {
      return m_semi_axes;
    }

    /// \brief Update the ellipsoid semi-axes (must be strictly positive).
    void setSemiAxes(const Vector3 & semi_axes)
    {
      m_semi_axes = semi_axes;
      check_and_init_axes();
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      if (this->joint1_id > 0)
        cdata.oMc1 = data.oMi[this->joint1_id] * this->joint1_placement;
      else
        cdata.oMc1 = this->joint1_placement;

      if (this->joint2_id > 0)
        cdata.oMc2 = data.oMi[this->joint2_id] * this->joint2_placement;
      else
        cdata.oMc2 = this->joint2_placement;

      cdata.d_world = cdata.oMc1.translation() - cdata.oMc2.translation();
      cdata.d_local.noalias() = cdata.oMc2.rotation().transpose() * cdata.d_world;

      cdata.D_d_local = m_inv_squared_axes.cwiseProduct(cdata.d_local);
      cdata.u_world.noalias() = cdata.oMc2.rotation() * cdata.D_d_local;

      cdata.constraint_residual[0] = cdata.d_local.dot(cdata.D_d_local) - Scalar(1);
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename JacobianMatrix>
    void jacobianImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JacobianMatrix> & _jacobian_matrix) const
    {
      JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();
      PINOCCHIO_CHECK_ARGUMENT_SIZE(jacobian_matrix.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(jacobian_matrix.cols(), model.nv);

      jacobian_matrix.setZero();

      const auto & sparsity1 = model.sparsity_pattern_vector[this->joint1_id];
      const auto & sparsity2 = model.sparsity_pattern_vector[this->joint2_id];

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const Vector3 & u_world = cdata.u_world;
      const Vector3 u_cross_d = u_world.cross(cdata.d_world);

      typedef DataTpl<Scalar, OtherOptions, JointCollectionTpl> DataType;

      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        const bool s1 = sparsity1[jj];
        const bool s2 = sparsity2[jj];
        if (!s1 && !s2)
          continue;

        typedef typename DataType::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);
        const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

        Scalar entry(0);

        if (s1)
        {
          const auto Jcol_local1 = oMc1.actInv(Jcol_motion);
          const Vector3 v_p1_world = oMc1.rotation() * Jcol_local1.linear();
          entry += Scalar(2) * u_world.dot(v_p1_world);
        }

        if (s2)
        {
          const auto Jcol_local2 = oMc2.actInv(Jcol_motion);
          const Vector3 v_pe_world = oMc2.rotation() * Jcol_local2.linear();
          const Vector3 omega_e_world = oMc2.rotation() * Jcol_local2.angular();
          entry -= Scalar(2) * u_world.dot(v_pe_world);
          entry += Scalar(2) * u_cross_d.dot(omega_e_world);
        }

        jacobian_matrix(0, jj) = entry;
      }
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix>
    typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
    jacobianMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianMatrixProductReturnType<InputMatrix>::type
        ReturnType;
      ReturnType res(residualSize(), mat.cols());
      jacobianMatrixProductImpl(model, data, cdata, mat.derived(), res, SetTo());
      return res;
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const
    {
      PINOCCHIO_UNUSED_VARIABLE(aot);
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), residualSize());

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      const auto & sparsity1 = model.sparsity_pattern_vector[this->joint1_id];
      const auto & sparsity2 = model.sparsity_pattern_vector[this->joint2_id];

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const Vector3 & u_world = cdata.u_world;
      const Vector3 u_cross_d = u_world.cross(cdata.d_world);

      typedef DataTpl<Scalar, OtherOptions, JointCollectionTpl> DataType;

      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        const bool s1 = sparsity1[jj];
        const bool s2 = sparsity2[jj];
        if (!s1 && !s2)
          continue;

        typedef typename DataType::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);
        const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

        Scalar entry(0);

        if (s1)
        {
          const auto Jcol_local1 = oMc1.actInv(Jcol_motion);
          const Vector3 v_p1_world = oMc1.rotation() * Jcol_local1.linear();
          entry += Scalar(2) * u_world.dot(v_p1_world);
        }

        if (s2)
        {
          const auto Jcol_local2 = oMc2.actInv(Jcol_motion);
          const Vector3 v_pe_world = oMc2.rotation() * Jcol_local2.linear();
          const Vector3 omega_e_world = oMc2.rotation() * Jcol_local2.angular();
          entry -= Scalar(2) * u_world.dot(v_pe_world);
          entry += Scalar(2) * u_cross_d.dot(omega_e_world);
        }

        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(0).noalias() -= entry * mat.row(jj);
        else
          res.row(0).noalias() += entry * mat.row(jj);
      }
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix>
    typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::type
    jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat) const
    {
      typedef typename traits<Self>::template JacobianTransposeMatrixProductReturnType<InputMatrix>::
        type ReturnType;
      ReturnType res(model.nv, mat.cols());
      jacobianTransposeMatrixProductImpl(model, data, cdata, mat.derived(), res, SetTo());
      return res;
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename InputMatrix,
      typename OutputMatrix,
      AssignmentOperatorType op>
    void jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const
    {
      PINOCCHIO_UNUSED_VARIABLE(aot);
      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      const auto & sparsity1 = model.sparsity_pattern_vector[this->joint1_id];
      const auto & sparsity2 = model.sparsity_pattern_vector[this->joint2_id];

      const SE3 & oMc1 = cdata.oMc1;
      const SE3 & oMc2 = cdata.oMc2;
      const Vector3 & u_world = cdata.u_world;
      const Vector3 u_cross_d = u_world.cross(cdata.d_world);

      typedef DataTpl<Scalar, OtherOptions, JointCollectionTpl> DataType;

      for (Eigen::Index jj = 0; jj < model.nv; ++jj)
      {
        const bool s1 = sparsity1[jj];
        const bool s2 = sparsity2[jj];
        if (!s1 && !s2)
          continue;

        typedef typename DataType::Matrix6x::ConstColXpr ConstColXpr;
        const ConstColXpr Jcol = data.J.col(jj);
        const MotionRef<const ConstColXpr> Jcol_motion(Jcol);

        Scalar entry(0);

        if (s1)
        {
          const auto Jcol_local1 = oMc1.actInv(Jcol_motion);
          const Vector3 v_p1_world = oMc1.rotation() * Jcol_local1.linear();
          entry += Scalar(2) * u_world.dot(v_p1_world);
        }

        if (s2)
        {
          const auto Jcol_local2 = oMc2.actInv(Jcol_motion);
          const Vector3 v_pe_world = oMc2.rotation() * Jcol_local2.linear();
          const Vector3 omega_e_world = oMc2.rotation() * Jcol_local2.angular();
          entry -= Scalar(2) * u_world.dot(v_pe_world);
          entry += Scalar(2) * u_cross_d.dot(omega_e_world);
        }

        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(jj).noalias() -= entry * mat.row(0);
        else
          res.row(jj).noalias() += entry * mat.row(0);
      }
    }

    template<
      int OtherOptions,
      int ForceOptions,
      template<typename, int> class JointCollectionTpl,
      typename ForceLike,
      typename ForceAllocator,
      ReferenceFrame rf>
    void mapConstraintForceToJointForcesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ForceLike> & constraint_forces,
      std::vector<ForceTpl<Scalar, ForceOptions>, ForceAllocator> & joint_forces,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_forces.size(), size_t(model.njoints));

      typedef ForceTpl<Scalar, ForceOptions> Force;

      const Scalar lambda = constraint_forces[0];
      const Vector3 F = Scalar(2) * lambda * cdata.u_world;

      if (this->joint1_id > 0)
      {
        const Vector3 r1 = cdata.oMc1.translation() - data.oMi[this->joint1_id].translation();
        Force f1(F, r1.cross(F));
        if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        {
          const Matrix3 & R1 = data.oMi[this->joint1_id].rotation();
          f1.linear() = R1.transpose() * f1.linear();
          f1.angular() = R1.transpose() * f1.angular();
        }
        joint_forces[this->joint1_id] += f1;
      }

      if (this->joint2_id > 0)
      {
        const Vector3 r12 = cdata.oMc1.translation() - data.oMi[this->joint2_id].translation();
        Force f2(-F, -r12.cross(F));
        if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        {
          const Matrix3 & R2 = data.oMi[this->joint2_id].rotation();
          f2.linear() = R2.transpose() * f2.linear();
          f2.angular() = R2.transpose() * f2.angular();
        }
        joint_forces[this->joint2_id] += f2;
      }
    }

    template<
      int OtherOptions,
      int MotionOptions,
      template<typename, int> class JointCollectionTpl,
      typename MotionAllocator,
      typename VectorLike,
      ReferenceFrame rf>
    void mapJointMotionsToConstraintMotionImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const std::vector<MotionTpl<Scalar, MotionOptions>, MotionAllocator> & joint_motions,
      const Eigen::MatrixBase<VectorLike> & constraint_motions_,
      ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions_.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.size(), size_t(model.njoints));

      auto & constraint_motions = constraint_motions_.const_cast_derived();

      const Vector3 & u_world = cdata.u_world;
      const Vector3 & d_world = cdata.d_world;
      const Vector3 u_cross_d = u_world.cross(d_world);

      Scalar value(0);

      if (this->joint1_id > 0)
      {
        Vector3 v_lin;
        Vector3 v_ang;
        if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        {
          const Matrix3 & R1 = data.oMi[this->joint1_id].rotation();
          v_lin.noalias() = R1 * joint_motions[this->joint1_id].linear();
          v_ang.noalias() = R1 * joint_motions[this->joint1_id].angular();
        }
        else
        {
          v_lin = joint_motions[this->joint1_id].linear();
          v_ang = joint_motions[this->joint1_id].angular();
        }
        const Vector3 r1 = cdata.oMc1.translation() - data.oMi[this->joint1_id].translation();
        const Vector3 v_p1_world = v_lin + v_ang.cross(r1);
        value += Scalar(2) * u_world.dot(v_p1_world);
      }

      if (this->joint2_id > 0)
      {
        Vector3 v_lin;
        Vector3 v_ang;
        if constexpr (std::is_same<ReferenceFrameTag<rf>, LocalFrameTag>::value)
        {
          const Matrix3 & R2 = data.oMi[this->joint2_id].rotation();
          v_lin.noalias() = R2 * joint_motions[this->joint2_id].linear();
          v_ang.noalias() = R2 * joint_motions[this->joint2_id].angular();
        }
        else
        {
          v_lin = joint_motions[this->joint2_id].linear();
          v_ang = joint_motions[this->joint2_id].angular();
        }
        const Vector3 r2 = cdata.oMc2.translation() - data.oMi[this->joint2_id].translation();
        const Vector3 v_pe_world = v_lin + v_ang.cross(r2);
        value -= Scalar(2) * u_world.dot(v_pe_world);
        value += Scalar(2) * u_cross_d.dot(v_ang);
      }

      constraint_motions[0] = value;
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename VectorNLike,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<VectorNLike> & diagonal_constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      PINOCCHIO_UNUSED_VARIABLE(diagonal_constraint_inertia);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "PointOnEllipsoidConstraintModel does not support appendCouplingConstraintInertias yet.");
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename MatrixOrMap,
      typename MapEnable,
      ReferenceFrame rf>
    void appendCouplingConstraintInertiasImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const internal::MatrixBlockElementTpl<MatrixOrMap, MapEnable> & constraint_inertia,
      const ReferenceFrameTag<rf> reference_frame) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);
      PINOCCHIO_UNUSED_VARIABLE(constraint_inertia);
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "PointOnEllipsoidConstraintModel does not support appendCouplingConstraintInertias yet.");
    }

  protected:
    void check_and_init_axes()
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        m_semi_axes[0] > Scalar(0) && m_semi_axes[1] > Scalar(0) && m_semi_axes[2] > Scalar(0),
        "PointOnEllipsoidConstraintModel: ellipsoid semi-axes must be strictly positive.");
      m_inv_squared_axes = m_semi_axes.cwiseProduct(m_semi_axes).cwiseInverse();
    }

    Vector3 m_semi_axes;
    Vector3 m_inv_squared_axes;
  };

  template<typename _Scalar, int _Options>
  struct PointOnEllipsoidConstraintDataTpl
  : ConstraintDataBase<PointOnEllipsoidConstraintDataTpl<_Scalar, _Options>>
  {
    typedef PointOnEllipsoidConstraintDataTpl Self;
    typedef ConstraintDataBase<Self> Base;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    typedef typename traits<Self>::ResidualVectorType ResidualVectorType;

    typedef SE3Tpl<Scalar, Options> SE3;
    typedef Eigen::Matrix<Scalar, 3, 1, Options> Vector3;

    using Base::classname;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    PointOnEllipsoidConstraintDataTpl()
    : constraint_residual(ResidualVectorType::Zero())
    , oMc1(SE3::Identity())
    , oMc2(SE3::Identity())
    , d_world(Vector3::Zero())
    , d_local(Vector3::Zero())
    , D_d_local(Vector3::Zero())
    , u_world(Vector3::Zero())
    {
    }

    explicit PointOnEllipsoidConstraintDataTpl(const ConstraintModel & cmodel)
    : constraint_residual(ResidualVectorType::Zero())
    , oMc1(SE3::Identity())
    , oMc2(SE3::Identity())
    , d_world(Vector3::Zero())
    , d_local(Vector3::Zero())
    , D_d_local(Vector3::Zero())
    , u_world(Vector3::Zero())
    {
      PINOCCHIO_UNUSED_VARIABLE(cmodel);
    }

    bool operator==(const PointOnEllipsoidConstraintDataTpl & other) const
    {
      return base() == other.base() && constraint_residual == other.constraint_residual
             && oMc1 == other.oMc1 && oMc2 == other.oMc2 && d_world == other.d_world
             && d_local == other.d_local && D_d_local == other.D_d_local
             && u_world == other.u_world;
    }

    bool operator!=(const PointOnEllipsoidConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    static std::string classnameImpl()
    {
      return std::string("PointOnEllipsoidConstraintData");
    }

    std::string shortnameImpl() const
    {
      return classname();
    }

    ResidualVectorType constraint_residual;
    SE3 oMc1;
    SE3 oMc2;
    Vector3 d_world;
    Vector3 d_local;
    Vector3 D_d_local;
    Vector3 u_world;
  };

} // namespace pinocchio
