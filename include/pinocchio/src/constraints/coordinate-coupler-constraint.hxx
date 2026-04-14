//
// Copyright (c) 2026 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/constraints.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/constraints.hpp"
#endif // PINOCCHIO_LSP

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <tuple>

namespace pinocchio
{

  // --------------------------------------------------------------
  // Cast
  // --------------------------------------------------------------
  template<typename NewScalar, typename Scalar, int Options>
  struct CastType<NewScalar, CoordinateCouplerConstraintModelTpl<Scalar, Options>>
  {
    typedef CoordinateCouplerConstraintModelTpl<NewScalar, Options> type;
  };

  // --------------------------------------------------------------
  // Traits
  // --------------------------------------------------------------
  template<typename _Scalar, int _Options>
  struct traits<CoordinateCouplerConstraintModelTpl<_Scalar, _Options>>
  {
    typedef CoordinateCouplerConstraintModelTpl<_Scalar, _Options> ConstraintModel;
    typedef CoordinateCouplerConstraintDataTpl<_Scalar, _Options> ConstraintData;

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
  struct traits<CoordinateCouplerConstraintDataTpl<_Scalar, _Options>>
  : traits<CoordinateCouplerConstraintModelTpl<_Scalar, _Options>>
  {
  };

  template<typename _Scalar, int _Options>
  struct CoordinateCouplerConstraintModelTpl
  : JointWiseConstraintModelBase<CoordinateCouplerConstraintModelTpl<_Scalar, _Options>>
  , ConstraintModelCommonParameters<CoordinateCouplerConstraintModelTpl<_Scalar, _Options>>
  {
    typedef CoordinateCouplerConstraintModelTpl Self;
    typedef JointWiseConstraintModelBase<Self> Base;
    typedef ConstraintModelCommonParameters<Self> BaseCommonParameters;
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
    static constexpr int SymmetricConeSize = traits<Self>::SymmetricConeSize;
    static constexpr int SymmetricConeScalingSize = traits<Self>::SymmetricConeScalingSize;

    typedef typename traits<Self>::ResidualVectorType ResidualVectorType;
    typedef typename traits<Self>::JacobianMatrixType JacobianMatrixType;
    typedef typename traits<Self>::ConeVectorType ConeVectorType;
    typedef typename traits<Self>::ConeScalingVectorType ConeScalingVectorType;
    typedef typename traits<Self>::RowVectorXs RowVectorXs;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      CompactTangentMap;
    typedef std::vector<JointIndex> JointIndexVector;

    template<typename NewScalar, int NewOptions>
    friend struct CoordinateCouplerConstraintModelTpl;

    template<typename NewScalar, int NewOptions>
    friend struct CoordinateCouplerConstraintDataTpl;

    using RootBase::classname;
    using RootBase::jacobianMatrixProduct;
    using RootBase::jacobianTransposeMatrixProduct;
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

    CoordinateCouplerConstraintModelTpl()
    : m_coordinate1_id(-1)
    , m_coordinate2_id(-1)
    , m_joint1_id(0)
    , m_joint2_id(0)
    , m_joint1_local_coordinate_id(-1)
    , m_joint2_local_coordinate_id(-1)
    , m_joint1_idx_v(0)
    , m_joint2_idx_v(0)
    , m_joint1_nq(0)
    , m_joint2_nq(0)
    , m_joint1_nv(0)
    , m_joint2_nv(0)
    , m_coordinate1_row_in_compact_tangent_map(-1)
    , m_coordinate2_row_in_compact_tangent_map(-1)
    , m_total_selected_nq(0)
    , m_max_selected_nv(0)
    {
      m_compliance = ResidualVectorType::Zero();
      m_baumgarte_parameters = BaumgarteCorrectorParameters();
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    explicit CoordinateCouplerConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    : CoordinateCouplerConstraintModelTpl(model, Eigen::Index(0), Eigen::Index(1))
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        model.nq >= 2, "CoordinateCouplerConstraintModel requires at least two coordinates.");
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    CoordinateCouplerConstraintModelTpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const Eigen::Index coordinate1_id,
      const Eigen::Index coordinate2_id)
    : Base(model)
    , m_coordinate1_id(coordinate1_id)
    , m_coordinate2_id(coordinate2_id)
    , m_joint1_id(0)
    , m_joint2_id(0)
    , m_joint1_local_coordinate_id(-1)
    , m_joint2_local_coordinate_id(-1)
    , m_joint1_idx_v(0)
    , m_joint2_idx_v(0)
    , m_joint1_nq(0)
    , m_joint2_nq(0)
    , m_joint1_nv(0)
    , m_joint2_nv(0)
    , m_coordinate1_row_in_compact_tangent_map(-1)
    , m_coordinate2_row_in_compact_tangent_map(-1)
    , m_total_selected_nq(0)
    , m_max_selected_nv(0)
    {
      init(model);
    }

    template<typename NewScalar>
    typename CastType<NewScalar, CoordinateCouplerConstraintModelTpl>::type cast() const
    {
      typedef typename CastType<NewScalar, CoordinateCouplerConstraintModelTpl>::type ReturnType;
      ReturnType res;
      Base::cast(res);
      BaseCommonParameters::template cast<NewScalar>(res);

      res.m_coordinate1_id = m_coordinate1_id;
      res.m_coordinate2_id = m_coordinate2_id;
      res.m_joint1_id = m_joint1_id;
      res.m_joint2_id = m_joint2_id;
      res.m_joint1_local_coordinate_id = m_joint1_local_coordinate_id;
      res.m_joint2_local_coordinate_id = m_joint2_local_coordinate_id;
      res.m_joint1_idx_v = m_joint1_idx_v;
      res.m_joint2_idx_v = m_joint2_idx_v;
      res.m_joint1_nq = m_joint1_nq;
      res.m_joint2_nq = m_joint2_nq;
      res.m_joint1_nv = m_joint1_nv;
      res.m_joint2_nv = m_joint2_nv;
      res.m_unique_joint_ids = m_unique_joint_ids;
      res.m_coordinate1_row_in_compact_tangent_map = m_coordinate1_row_in_compact_tangent_map;
      res.m_coordinate2_row_in_compact_tangent_map = m_coordinate2_row_in_compact_tangent_map;
      res.m_total_selected_nq = m_total_selected_nq;
      res.m_max_selected_nv = m_max_selected_nv;

      return res;
    }

    bool operator==(const CoordinateCouplerConstraintModelTpl & other) const
    {
      return base() == other.base() && base_common_parameters() == other.base_common_parameters()
             && m_coordinate1_id == other.m_coordinate1_id
             && m_coordinate2_id == other.m_coordinate2_id && m_joint1_id == other.m_joint1_id
             && m_joint2_id == other.m_joint2_id
             && m_joint1_local_coordinate_id == other.m_joint1_local_coordinate_id
             && m_joint2_local_coordinate_id == other.m_joint2_local_coordinate_id
             && m_joint1_idx_v == other.m_joint1_idx_v && m_joint2_idx_v == other.m_joint2_idx_v
             && m_joint1_nq == other.m_joint1_nq && m_joint2_nq == other.m_joint2_nq
             && m_joint1_nv == other.m_joint1_nv && m_joint2_nv == other.m_joint2_nv
             && m_unique_joint_ids == other.m_unique_joint_ids
             && m_coordinate1_row_in_compact_tangent_map
                  == other.m_coordinate1_row_in_compact_tangent_map
             && m_coordinate2_row_in_compact_tangent_map
                  == other.m_coordinate2_row_in_compact_tangent_map
             && m_total_selected_nq == other.m_total_selected_nq
             && m_max_selected_nv == other.m_max_selected_nv;
    }

    bool operator!=(const CoordinateCouplerConstraintModelTpl & other) const
    {
      return !(*this == other);
    }

    static std::string classnameImpl()
    {
      return std::string("CoordinateCouplerConstraintModel");
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

    Eigen::Index getCoordinate1Id() const
    {
      return m_coordinate1_id;
    }

    Eigen::Index getCoordinate2Id() const
    {
      return m_coordinate2_id;
    }

    JointIndex getJoint1Id() const
    {
      return m_joint1_id;
    }

    JointIndex getJoint2Id() const
    {
      return m_joint2_id;
    }

    int getJoint1Nv() const
    {
      return m_joint1_nv;
    }

    int getJoint2Nv() const
    {
      return m_joint2_nv;
    }

    bool couplesDifferentJoints() const
    {
      return m_joint1_id != m_joint2_id;
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void calcImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      ConstraintData & cdata) const
    {
      PINOCCHIO_UNUSED_VARIABLE(model);

      cdata.compact_tangent_map.resize(m_total_selected_nq, m_max_selected_nv);
      cdata.coordinate1_tangent_map.resize(m_joint1_nv);
      cdata.coordinate2_tangent_map.resize(m_joint2_nv);

      pinocchio::compactTangentMap(model, m_unique_joint_ids, data.q_in, cdata.compact_tangent_map);

      cdata.constraint_residual[0] = data.q_in[m_coordinate1_id] - data.q_in[m_coordinate2_id];

      cdata.coordinate1_tangent_map =
        cdata.compact_tangent_map.row(m_coordinate1_row_in_compact_tangent_map).head(m_joint1_nv);
      cdata.coordinate2_tangent_map =
        cdata.compact_tangent_map.row(m_coordinate2_row_in_compact_tangent_map).head(m_joint2_nv);
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void getRowSparsityPatternImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id,
      BooleanVector & result) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id == 0);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);

      result = model.sparsity_pattern_vector[m_joint1_id];
      if (m_joint2_id != m_joint1_id)
        result = result.array() || model.sparsity_pattern_vector[m_joint2_id].array();
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void getRowIndexesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::Index row_id,
      EigenIndexVector & result) const
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(row_id == 0);
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(cdata);

      if (m_joint1_id == m_joint2_id)
      {
        result = model.span_indexes_vector[m_joint1_id];
        return;
      }

      const auto & indexes1 = model.span_indexes_vector[m_joint1_id];
      const auto & indexes2 = model.span_indexes_vector[m_joint2_id];

      result.clear();
      result.reserve(indexes1.size() + indexes2.size());
      std::set_union(
        indexes1.begin(), indexes1.end(), indexes2.begin(), indexes2.end(),
        std::back_inserter(result));
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
      PINOCCHIO_UNUSED_VARIABLE(data);

      JacobianMatrix & jacobian_matrix = _jacobian_matrix.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(jacobian_matrix.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(jacobian_matrix.cols(), model.nv);

      jacobian_matrix.setZero();

      if (m_joint1_id == m_joint2_id)
      {
        jacobian_matrix.row(0).segment(m_joint1_idx_v, m_joint1_nv) =
          cdata.coordinate1_tangent_map - cdata.coordinate2_tangent_map;
      }
      else
      {
        jacobian_matrix.row(0).segment(m_joint1_idx_v, m_joint1_nv) = cdata.coordinate1_tangent_map;
        jacobian_matrix.row(0).segment(m_joint2_idx_v, m_joint2_nv) =
          -cdata.coordinate2_tangent_map;
      }
    }

    template<
      int OtherOptions,
      typename InputMatrix,
      template<typename, int> class JointCollectionTpl>
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
      jacobianMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    template<
      int OtherOptions,
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op>
    void jacobianMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const
    {
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(aot);

      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), model.nv);
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.cols(), res.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), residualSize());

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      const auto term1 = cdata.coordinate1_tangent_map * mat.middleRows(m_joint1_idx_v, m_joint1_nv);

      if (m_joint1_id == m_joint2_id)
      {
        const auto term2 =
          cdata.coordinate2_tangent_map * mat.middleRows(m_joint1_idx_v, m_joint1_nv);
        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(0).noalias() -= term1 - term2;
        else
          res.row(0).noalias() += term1 - term2;
      }
      else
      {
        const auto term2 =
          cdata.coordinate2_tangent_map * mat.middleRows(m_joint2_idx_v, m_joint2_nv);
        if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
          res.row(0).noalias() -= term1 - term2;
        else
          res.row(0).noalias() += term1 - term2;
      }
    }

    template<
      int OtherOptions,
      typename InputMatrix,
      template<typename, int> class JointCollectionTpl>
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
      jacobianTransposeMatrixProduct(model, data, cdata, mat.derived(), res);
      return res;
    }

    template<
      int OtherOptions,
      typename InputMatrix,
      typename OutputMatrix,
      template<typename, int> class JointCollectionTpl,
      AssignmentOperatorType op>
    void jacobianTransposeMatrixProductImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<InputMatrix> & mat,
      const Eigen::MatrixBase<OutputMatrix> & _res,
      AssignmentOperatorTag<op> aot) const
    {
      PINOCCHIO_UNUSED_VARIABLE(data);
      PINOCCHIO_UNUSED_VARIABLE(aot);

      OutputMatrix & res = _res.const_cast_derived();

      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), mat.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), model.nv);

      if constexpr (std::is_same<AssignmentOperatorTag<op>, SetTo>::value)
        res.setZero();

      if constexpr (std::is_same<AssignmentOperatorTag<op>, RmTo>::value)
      {
        res.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() -=
          cdata.coordinate1_tangent_map.transpose() * mat.row(0);
        if (m_joint1_id == m_joint2_id)
        {
          res.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() +=
            cdata.coordinate2_tangent_map.transpose() * mat.row(0);
        }
        else
        {
          res.middleRows(m_joint2_idx_v, m_joint2_nv).noalias() +=
            cdata.coordinate2_tangent_map.transpose() * mat.row(0);
        }
      }
      else
      {
        res.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() +=
          cdata.coordinate1_tangent_map.transpose() * mat.row(0);
        if (m_joint1_id == m_joint2_id)
        {
          res.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() -=
            cdata.coordinate2_tangent_map.transpose() * mat.row(0);
        }
        else
        {
          res.middleRows(m_joint2_idx_v, m_joint2_nv).noalias() -=
            cdata.coordinate2_tangent_map.transpose() * mat.row(0);
        }
      }
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename ConstraintForcesLike,
      typename JointTorquesLike>
    void mapConstraintForceToJointTorquesImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<ConstraintForcesLike> & constraint_forces,
      const Eigen::MatrixBase<JointTorquesLike> & joint_torques_) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_forces.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_torques_.rows(), model.nv);
      PINOCCHIO_UNUSED_VARIABLE(data);

      auto & joint_torques = joint_torques_.const_cast_derived();

      joint_torques.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() +=
        cdata.coordinate1_tangent_map.transpose() * constraint_forces.row(0);

      if (m_joint1_id == m_joint2_id)
      {
        joint_torques.middleRows(m_joint1_idx_v, m_joint1_nv).noalias() -=
          cdata.coordinate2_tangent_map.transpose() * constraint_forces.row(0);
      }
      else
      {
        joint_torques.middleRows(m_joint2_idx_v, m_joint2_nv).noalias() -=
          cdata.coordinate2_tangent_map.transpose() * constraint_forces.row(0);
      }
    }

    template<
      int OtherOptions,
      template<typename, int> class JointCollectionTpl,
      typename JointMotionsLike,
      typename ConstraintMotionsLike>
    void mapJointMotionsToConstraintMotionImpl(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const DataTpl<Scalar, OtherOptions, JointCollectionTpl> & data,
      const ConstraintData & cdata,
      const Eigen::MatrixBase<JointMotionsLike> & joint_motions,
      const Eigen::MatrixBase<ConstraintMotionsLike> & constraint_motions_) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(constraint_motions_.rows(), residualSize());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(joint_motions.rows(), model.nv);
      PINOCCHIO_UNUSED_VARIABLE(data);

      auto & constraint_motions = constraint_motions_.const_cast_derived();
      constraint_motions.row(0).noalias() =
        cdata.coordinate1_tangent_map * joint_motions.middleRows(m_joint1_idx_v, m_joint1_nv);

      if (m_joint1_id == m_joint2_id)
      {
        constraint_motions.row(0).noalias() -=
          cdata.coordinate2_tangent_map * joint_motions.middleRows(m_joint1_idx_v, m_joint1_nv);
      }
      else
      {
        constraint_motions.row(0).noalias() -=
          cdata.coordinate2_tangent_map * joint_motions.middleRows(m_joint2_idx_v, m_joint2_nv);
      }
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
      PINOCCHIO_UNUSED_VARIABLE(reference_frame);

      PINOCCHIO_CHECK_ARGUMENT_SIZE(diagonal_constraint_inertia.size(), residualSize());

      if (m_joint1_id != m_joint2_id)
      {
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument,
          "CoordinateCouplerConstraintModel only supports appendCouplingConstraintInertias when "
          "both coordinates belong to the same joint.");
      }

      const Scalar inertia_value = diagonal_constraint_inertia.derived().coeff(0);
      const RowVectorXs tangent =
        cdata.coordinate1_tangent_map - cdata.coordinate2_tangent_map;

      data.joint_apparent_inertia[m_joint1_id].noalias() +=
        inertia_value * tangent.transpose() * tangent;
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
      assert(constraint_inertia.size() == residualSize());

      switch (constraint_inertia.type())
      {
      case internal::MatrixBlockType::Zero: {
        break;
      }
      case internal::MatrixBlockType::Identity: {
        appendCouplingConstraintInertiasImpl(
          model, data, cdata, ResidualVectorType::Ones(), reference_frame);
        break;
      }
      case internal::MatrixBlockType::ScalarIdentity: {
        appendCouplingConstraintInertiasImpl(
          model, data, cdata,
          ResidualVectorType::Constant(constraint_inertia.container()(0, 0)), reference_frame);
        break;
      }
      case internal::MatrixBlockType::Diagonal: {
        ResidualVectorType diagonal;
        diagonal[0] = constraint_inertia.container()(0, 0);
        appendCouplingConstraintInertiasImpl(model, data, cdata, diagonal, reference_frame);
        break;
      }
      case internal::MatrixBlockType::Plain: {
        ResidualVectorType inertia_value;
        inertia_value[0] = constraint_inertia.container()(0, 0);
        appendCouplingConstraintInertiasImpl(model, data, cdata, inertia_value, reference_frame);
        break;
      }
      default:
        assert(false && "Invalid MatrixBlockType for CoordinateCouplerConstraintModel.");
        PINOCCHIO_THROW_PRETTY(
          std::invalid_argument, "Invalid MatrixBlockType for CoordinateCouplerConstraintModel.");
      }
    }

  protected:
    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    void init(const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model)
    {
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        m_coordinate1_id >= 0 && m_coordinate1_id < model.nq,
        "coordinate1_id should belong to [0, model.nq).");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        m_coordinate2_id >= 0 && m_coordinate2_id < model.nq,
        "coordinate2_id should belong to [0, model.nq).");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(
        m_coordinate1_id != m_coordinate2_id,
        "CoordinateCouplerConstraintModel requires two distinct coordinates.");

      std::tie(m_joint1_id, m_joint1_local_coordinate_id, m_joint1_nq, m_joint1_nv, m_joint1_idx_v)
        = findJointForCoordinate(model, m_coordinate1_id);
      std::tie(m_joint2_id, m_joint2_local_coordinate_id, m_joint2_nq, m_joint2_nv, m_joint2_idx_v)
        = findJointForCoordinate(model, m_coordinate2_id);

      PINOCCHIO_CHECK_INPUT_ARGUMENT(m_joint1_nv > 0, "coordinate1_id belongs to a zero-DoF joint.");
      PINOCCHIO_CHECK_INPUT_ARGUMENT(m_joint2_nv > 0, "coordinate2_id belongs to a zero-DoF joint.");

      m_unique_joint_ids.clear();
      m_unique_joint_ids.push_back(m_joint1_id);
      m_coordinate1_row_in_compact_tangent_map = m_joint1_local_coordinate_id;

      if (m_joint2_id == m_joint1_id)
      {
        m_coordinate2_row_in_compact_tangent_map = m_joint2_local_coordinate_id;
        m_total_selected_nq = m_joint1_nq;
        m_max_selected_nv = m_joint1_nv;
      }
      else
      {
        m_unique_joint_ids.push_back(m_joint2_id);
        m_coordinate2_row_in_compact_tangent_map = m_joint1_nq + m_joint2_local_coordinate_id;
        m_total_selected_nq = m_joint1_nq + m_joint2_nq;
        m_max_selected_nv = std::max(m_joint1_nv, m_joint2_nv);
      }

      m_compliance = ResidualVectorType::Zero();
      m_baumgarte_parameters = BaumgarteCorrectorParameters();
    }

    template<int OtherOptions, template<typename, int> class JointCollectionTpl>
    static std::tuple<JointIndex, Eigen::Index, int, int, int> findJointForCoordinate(
      const ModelTpl<Scalar, OtherOptions, JointCollectionTpl> & model,
      const Eigen::Index coordinate_id)
    {
      for (JointIndex joint_id = 1; joint_id < JointIndex(model.njoints); ++joint_id)
      {
        const Eigen::Index idx_q = model.idx_qs[joint_id];
        const Eigen::Index nq = model.nqs[joint_id];
        if (coordinate_id >= idx_q && coordinate_id < idx_q + nq)
        {
          return std::make_tuple(
            joint_id, coordinate_id - idx_q, model.nqs[joint_id], model.nvs[joint_id],
            model.idx_vs[joint_id]);
        }
      }

      PINOCCHIO_THROW_PRETTY(
        std::invalid_argument,
        "Could not find a joint associated to the input configuration coordinate.");
      PINOCCHIO_UNREACHABLE();
    }

    using BaseCommonParameters::m_baumgarte_parameters;
    using BaseCommonParameters::m_compliance;

    Eigen::Index m_coordinate1_id;
    Eigen::Index m_coordinate2_id;

    JointIndex m_joint1_id;
    JointIndex m_joint2_id;

    Eigen::Index m_joint1_local_coordinate_id;
    Eigen::Index m_joint2_local_coordinate_id;

    int m_joint1_idx_v;
    int m_joint2_idx_v;
    int m_joint1_nq;
    int m_joint2_nq;
    int m_joint1_nv;
    int m_joint2_nv;

    JointIndexVector m_unique_joint_ids;
    Eigen::Index m_coordinate1_row_in_compact_tangent_map;
    Eigen::Index m_coordinate2_row_in_compact_tangent_map;
    int m_total_selected_nq;
    int m_max_selected_nv;
  };

  template<typename _Scalar, int _Options>
  struct CoordinateCouplerConstraintDataTpl
  : ConstraintDataBase<CoordinateCouplerConstraintDataTpl<_Scalar, _Options>>
  {
    typedef CoordinateCouplerConstraintDataTpl Self;
    typedef ConstraintDataBase<Self> Base;

    typedef typename traits<Self>::ConstraintModel ConstraintModel;
    typedef typename traits<Self>::ConstraintData ConstraintData;

    typedef typename traits<Self>::Scalar Scalar;
    static constexpr int Options = traits<Self>::Options;

    typedef typename traits<Self>::ResidualVectorType ResidualVectorType;
    typedef typename traits<Self>::RowVectorXs RowVectorXs;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      CompactTangentMap;

    using Base::classname;

    Base & base()
    {
      return static_cast<Base &>(*this);
    }

    const Base & base() const
    {
      return static_cast<const Base &>(*this);
    }

    CoordinateCouplerConstraintDataTpl()
    : constraint_residual(ResidualVectorType::Zero())
    {
    }

    explicit CoordinateCouplerConstraintDataTpl(const ConstraintModel & cmodel)
    : constraint_residual(ResidualVectorType::Zero())
    {
      compact_tangent_map.resize(cmodel.m_total_selected_nq, cmodel.m_max_selected_nv);
      coordinate1_tangent_map.resize(cmodel.m_joint1_nv);
      coordinate2_tangent_map.resize(cmodel.m_joint2_nv);
      compact_tangent_map.setZero();
      coordinate1_tangent_map.setZero();
      coordinate2_tangent_map.setZero();
    }

    bool operator==(const CoordinateCouplerConstraintDataTpl & other) const
    {
      return base() == other.base() && constraint_residual == other.constraint_residual
             && compact_tangent_map == other.compact_tangent_map
             && coordinate1_tangent_map == other.coordinate1_tangent_map
             && coordinate2_tangent_map == other.coordinate2_tangent_map;
    }

    bool operator!=(const CoordinateCouplerConstraintDataTpl & other) const
    {
      return !(*this == other);
    }

    static std::string classnameImpl()
    {
      return std::string("CoordinateCouplerConstraintData");
    }

    std::string shortnameImpl() const
    {
      return classname();
    }

    ResidualVectorType constraint_residual;
    CompactTangentMap compact_tangent_map;
    RowVectorXs coordinate1_tangent_map;
    RowVectorXs coordinate2_tangent_map;
  };

} // namespace pinocchio
