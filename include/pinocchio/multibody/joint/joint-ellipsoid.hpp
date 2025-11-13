//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_multibody_joint_ellipsoid_hpp__
#define __pinocchio_multibody_joint_ellipsoid_hpp__

#include "pinocchio/macros.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/joint-motion-subspace.hpp"
#include "pinocchio/math/sincos.hpp"
#include "pinocchio/math/matrix.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/multibody/joint-motion-subspace.hpp"

namespace pinocchio
{
  template<typename Scalar, int Options>
  struct JointEllipsoidTpl;

  template<typename _Scalar, int _Options>
  struct traits<JointEllipsoidTpl<_Scalar, _Options>>
  {
    enum
    {
      NQ = 3,
      NV = 3,
      NVExtended = 3
    };
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options
    };
    typedef JointDataEllipsoidTpl<Scalar, Options> JointDataDerived;
    typedef JointModelEllipsoidTpl<Scalar, Options> JointModelDerived;
    typedef JointMotionSubspaceTpl<3, Scalar, Options, 3> Constraint_t;
    typedef SE3Tpl<Scalar, Options> Transformation_t;

    typedef MotionTpl<Scalar, Options> Motion_t;
    typedef MotionTpl<Scalar, Options> Bias_t;

    // [ABA]
    typedef Eigen::Matrix<Scalar, 6, NV, Options> U_t;
    typedef Eigen::Matrix<Scalar, NV, NV, Options> D_t;
    typedef Eigen::Matrix<Scalar, 6, NV, Options> UD_t;

    typedef Eigen::Matrix<Scalar, NQ, 1, Options> ConfigVector_t;
    typedef Eigen::Matrix<Scalar, NV, 1, Options> TangentVector_t;

    typedef boost::mpl::false_ is_mimicable_t; // not mimicable

    PINOCCHIO_JOINT_DATA_BASE_ACCESSOR_DEFAULT_RETURN_TYPE
  };

  template<typename _Scalar, int _Options>
  struct traits<JointDataEllipsoidTpl<_Scalar, _Options>>
  {
    typedef JointEllipsoidTpl<_Scalar, _Options> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options>
  struct traits<JointModelEllipsoidTpl<_Scalar, _Options>>
  {
    typedef JointEllipsoidTpl<_Scalar, _Options> JointDerived;
    typedef _Scalar Scalar;
  };

  template<typename _Scalar, int _Options>
  struct JointDataEllipsoidTpl : public JointDataBase<JointDataEllipsoidTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef JointEllipsoidTpl<_Scalar, _Options> JointDerived;
    PINOCCHIO_JOINT_DATA_TYPEDEF_TEMPLATE(JointDerived);
    PINOCCHIO_JOINT_DATA_BASE_DEFAULT_ACCESSOR

    ConfigVector_t joint_q;
    TangentVector_t joint_v;

    Constraint_t S;
    Constraint_t Sdot;
    Transformation_t M;
    Motion_t v;
    Bias_t c;

    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;
    D_t StU;

    JointDataEllipsoidTpl()
    : joint_q(ConfigVector_t::Zero())
    , joint_v(TangentVector_t::Zero())
    , S()
    , Sdot()
    , M(Transformation_t::Identity())
    , v(Motion_t::Zero())
    , c(Bias_t::Zero())
    , U(U_t::Zero())
    , Dinv(D_t::Zero())
    , UDinv(UD_t::Zero())
    , StU(D_t::Zero())
    {
    }

    static std::string classname()
    {
      return std::string("JointDataEllipsoid");
    }
    std::string shortname() const
    {
      return classname();
    }

  }; // struct JointDataEllipsoidTpl

  PINOCCHIO_JOINT_CAST_TYPE_SPECIALIZATION(JointModelEllipsoidTpl);
  template<typename _Scalar, int _Options>
  struct JointModelEllipsoidTpl : public JointModelBase<JointModelEllipsoidTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef JointEllipsoidTpl<_Scalar, _Options> JointDerived;
    PINOCCHIO_JOINT_TYPEDEF_TEMPLATE(JointDerived);

    typedef JointModelBase<JointModelEllipsoidTpl> Base;
    using Base::id;
    using Base::idx_q;
    using Base::idx_v;
    using Base::idx_vExtended;
    using Base::setIndexes;

    typedef typename Transformation_t::Vector3 Vector3;

    Scalar radius_a;
    Scalar radius_b;
    Scalar radius_c;

    JointDataDerived createData() const
    {
      return JointDataDerived();
    }

    JointModelEllipsoidTpl()
    : radius_a(Scalar(0.01))
    , radius_b(Scalar(0.01))
    , radius_c(Scalar(0.01))
    {
    }

    JointModelEllipsoidTpl(const Scalar & a, const Scalar & b, const Scalar & c)
    : radius_a(a)
    , radius_b(b)
    , radius_c(c)
    {
    }

    const std::vector<bool> hasConfigurationLimit() const
    {
      return {true, true, true};
    }

    const std::vector<bool> hasConfigurationLimitInTangent() const
    {
      return {true, true, true};
    }

    void computeSpatialTransform(
      Scalar c0, Scalar s0, Scalar c1, Scalar s1, Scalar c2, Scalar s2, JointDataDerived & data) const
    {
      // clang-format off
      data.M.rotation() << c1 * c2, -c1 * s2, s1, c0 * s2 + c2 * s0 * s1, c0 * c2 - s0 * s1 * s2,
      -c1 * s0, -c0 * c2 * s1 + s0 * s2, c0 * s1 * s2 + c2 * s0, c0 * c1;
      // clang-format on
      Scalar nx, ny, nz;
      nx = s1;
      ny = -s0 * c1;
      nz = c0 * c1;

      data.M.translation() << radius_a * nx, radius_b * ny, radius_c * nz;

    }

    template<typename ConfigVector>
    void calc(JointDataDerived & data, const typename Eigen::MatrixBase<ConfigVector> & qs) const
    {
      data.joint_q = qs.template segment<NQ>(idx_q());

      Scalar c0, s0;
      SINCOS(data.joint_q(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(data.joint_q(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(data.joint_q(2), &s2, &c2);

      computeSpatialTransform(c0,s0,c1,s1,c2,s2, data);

      computeMotionSubspace(s0, c0, s1, c1, s2, c2, data);
    }

    template<typename TangentVector>
    void
    calc(JointDataDerived & data, const Blank, const typename Eigen::MatrixBase<TangentVector> & vs)
      const
    {
      data.joint_v = vs.template segment<NV>(idx_v());

      // Compute S from already-set joint_q
      Scalar c0, s0;
      SINCOS(data.joint_q(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(data.joint_q(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(data.joint_q(2), &s2, &c2);

      Scalar dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0, dndotz_dqdot1;
      dndotx_dqdot1 = c1;
      dndoty_dqdot0 = -c0 * c1;
      dndoty_dqdot1 = s0 * s1;
      dndotz_dqdot0 = -c1 * s0;
      dndotz_dqdot1 = -c0 * s1;

      computeMotionSubspace(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);

      // Velocity part
      data.v.toVector().noalias() = data.S.matrix() * data.joint_v;

      computeMotionSubspaceDerivative(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);

      // data.c.toVector().noalias() = data.Sdot.matrix() * data.joint_v;
      computeBiais(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);
    }

    template<typename ConfigVector, typename TangentVector>
    void calc(
      JointDataDerived & data,
      const typename Eigen::MatrixBase<ConfigVector> & qs,
      const typename Eigen::MatrixBase<TangentVector> & vs) const
    {
      data.joint_q = qs.template segment<NQ>(idx_q());

      Scalar c0, s0;
      SINCOS(data.joint_q(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(data.joint_q(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(data.joint_q(2), &s2, &c2);

      computeSpatialTransform(c0, s0, c1, s1, c2, s2, data);
      
      Scalar dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0, dndotz_dqdot1;
      dndotx_dqdot1 = c1;
      dndoty_dqdot0 = -c0 * c1;
      dndoty_dqdot1 = s0 * s1;
      dndotz_dqdot0 = -c1 * s0;
      dndotz_dqdot1 = -c0 * s1;

      computeMotionSubspace(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);

      // Velocity part
      data.joint_v = vs.template segment<NV>(idx_v());
      data.v.toVector().noalias() = data.S.matrix() * data.joint_v;

      computeMotionSubspaceDerivative(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);

      // data.c.toVector().noalias() = data.Sdot.matrix() * data.joint_v;
      computeBiais(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);
    }

    template<typename VectorLike, typename Matrix6Like>
    void calc_aba(
      JointDataDerived & data,
      const Eigen::MatrixBase<VectorLike> & armature,
      const Eigen::MatrixBase<Matrix6Like> & I,
      const bool update_I) const
    {
      data.U.noalias() = I * data.S.matrix();
      data.StU.noalias() = data.S.transpose() * data.U;
      data.StU.diagonal() += armature;
      internal::PerformStYSInversion<Scalar>::run(data.StU, data.Dinv);

      data.UDinv.noalias() = data.U * data.Dinv;

      if (update_I)
        PINOCCHIO_EIGEN_CONST_CAST(Matrix6Like, I).noalias() -= data.UDinv * data.U.transpose();
    }

    static std::string classname()
    {
      return std::string("JointModelEllipsoid");
    }
    std::string shortname() const
    {
      return classname();
    }

    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    JointModelEllipsoidTpl<NewScalar, Options> cast() const
    {
      typedef JointModelEllipsoidTpl<NewScalar, Options> ReturnType;
      ReturnType res;
      res.setIndexes(id(), idx_q(), idx_v(), idx_vExtended());
      return res;
    }

    template<typename ConfigVector>
    Vector3 computeTranslations(const Eigen::MatrixBase<ConfigVector> & qs) const
    {
      Scalar c0, s0;
      SINCOS(qs(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(qs(1), &s1, &c1);

      return computeTranslations(s0, c0, s1, c1);
    }

    Vector3 computeTranslations(
      const Scalar & s0, const Scalar & c0, const Scalar & s1, const Scalar & c1) const
    {
      Scalar nx, ny, nz;
      nx = s1;
      ny = -s0 * c1;
      nz = c0 * c1;

      return Vector3(radius_a * nx, radius_b * ny, radius_c * nz);
    }

    template<typename ConfigVector, typename TangentVector>
    Vector3 computeTranslationVelocities(
      const Eigen::MatrixBase<ConfigVector> & qs,
      const Eigen::MatrixBase<TangentVector> & vs) const
    {
      Scalar c0, s0;
      SINCOS(qs(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(qs(1), &s1, &c1);

      return computeTranslationVelocities(s0, c0, s1, c1, vs(0), vs(1));
    }

    Vector3 computeTranslationVelocities(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & q0dot,
      const Scalar & q1dot) const
    {
      Vector3 v;
      v(0) = radius_a * c1 * q1dot;
      v(1) = radius_b * (-c0 * c1 * q0dot + s0 * s1 * q1dot);
      v(2) = radius_c * (-s0 * c1 * q0dot - c0 * s1 * q1dot);
      return v;
    }

    template<typename ConfigVector, typename TangentVector, typename TangentVector2>
    Vector3 computeTranslationAccelerations(
      const Eigen::MatrixBase<ConfigVector> & qs,
      const Eigen::MatrixBase<TangentVector> & vs,
      const Eigen::MatrixBase<TangentVector2> & as) const
    {
      Vector3 a;
      Scalar c0, s0;
      SINCOS(qs(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(qs(1), &s1, &c1);
      return computeTranslationAccelerations(s0, c0, s1, c1, vs(0), vs(1), as(0), as(1));
    }

    Vector3 computeTranslationAccelerations(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & q0dot,
      const Scalar & q1dot,
      const Scalar & q0ddot,
      const Scalar & q1ddot) const
    {
      Vector3 a;
      a(0) = radius_a * (-s1 * q1dot * q1dot + c1 * q1ddot);
      a(1) =
        radius_b
        * (s0 * c1 * q0dot * q0dot + c0 * s1 * q0dot * q1dot - c0 * c1 * q0ddot + c0 * s1 * q1dot * q0dot + s0 * c1 * q1dot * q1dot + s0 * s1 * q1ddot);
      a(2) =
        radius_c
        * (-c0 * c1 * q0dot * q0dot + s0 * s1 * q0dot * q1dot - s0 * c1 * q0ddot + s0 * s1 * q1dot * q0dot - c0 * c1 * q1dot * q1dot - c0 * s1 * q1ddot);
      return a;
    }

    void computeMotionSubspace(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & s2,
      const Scalar & c2,
      JointDataDerived & data) const
    {

      Scalar dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0, dndotz_dqdot1;
      dndotx_dqdot1 = c1;
      dndoty_dqdot0 = -c0 * c1;
      dndoty_dqdot1 = s0 * s1;
      dndotz_dqdot0 = -c1 * s0;
      dndotz_dqdot1 = -c0 * s1;

      computeMotionSubspace(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);
    }
    void computeMotionSubspace(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & s2,
      const Scalar & c2,
      const Scalar & dndotx_dqdot1,
      const Scalar & dndoty_dqdot0,
      const Scalar & dndoty_dqdot1,
      const Scalar & dndotz_dqdot0,
      const Scalar & dndotz_dqdot1,
      JointDataDerived & data) const
    {
      Scalar S_11, S_21, S_31, S_12, S_22, S_32;
      // clang-format off
      S_11 = dndoty_dqdot0 * radius_b * (c0 * s2 + c2 * s0 * s1)
             + dndotz_dqdot0 * radius_c * (-c0 * c2 * s1 + s0 * s2);

      S_21 = -dndoty_dqdot0 * radius_b * (-c0 * c2 + s0 * s1 * s2)
             + dndotz_dqdot0 * radius_c * (c0 * s1 * s2 + c2 * s0);

      S_31 = c1 * (-dndoty_dqdot0 * radius_b * s0 + dndotz_dqdot0 * radius_c * c0);

      S_12 = dndotx_dqdot1 * radius_a * c1 * c2
             + dndoty_dqdot1 * radius_b * (c0 * s2 + c2 * s0 * s1)
             + dndotz_dqdot1 * radius_c * (-c0 * c2 * s1 + s0 * s2);

      S_22 = -dndotx_dqdot1 * radius_a * c1 * s2
             - dndoty_dqdot1 * radius_b * (-c0 * c2 + s0 * s1 * s2)
             + dndotz_dqdot1 * radius_c * (c0 * s1 * s2 + c2 * s0);

      S_32 = dndotx_dqdot1 * radius_a * s1 - dndoty_dqdot1 * radius_b * c1 * s0
             + dndotz_dqdot1 * radius_c * c0 * c1;

      data.S.matrix() << S_11, S_12, Scalar(0), S_21, S_22, Scalar(0), S_31, S_32, Scalar(0),
        c1 * c2, s2, Scalar(0), -c1 * s2, c2, Scalar(0), s1, Scalar(0), Scalar(1);
      // clang-format on
    }

    template<typename ConfigVector, typename TangentVector>
    void computeMotionSubspaceDerivative(
      JointDataDerived & data,
      const Eigen::MatrixBase<ConfigVector> & qs,
      const Eigen::MatrixBase<TangentVector> & vs) const
    {
      Scalar c0, s0;
      SINCOS(qs(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(qs(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(qs(2), &s2, &c2);

      Scalar dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0, dndotz_dqdot1;
      dndotx_dqdot1 = c1;
      dndoty_dqdot0 = -c0 * c1;
      dndoty_dqdot1 = s0 * s1;
      dndotz_dqdot0 = -c1 * s0;
      dndotz_dqdot1 = -c0 * s1;

      computeMotionSubspaceDerivative(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);
    }

    void computeMotionSubspaceDerivative(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & s2,
      const Scalar & c2,
      const Scalar & dndotx_dqdot1,
      const Scalar & dndoty_dqdot0,
      const Scalar & dndoty_dqdot1,
      const Scalar & dndotz_dqdot0,
      const Scalar & dndotz_dqdot1,
      JointDataDerived & data) const
    {
      // Compute Sdot for bias acceleration
      Scalar qdot0, qdot1, qdot2;
      qdot0 = data.joint_v(0);
      qdot1 = data.joint_v(1);
      qdot2 = data.joint_v(2);

      Scalar Sdot_11, Sdot_21, Sdot_31, Sdot_41, Sdot_51, Sdot_61;
      Scalar Sdot_12, Sdot_22, Sdot_32, Sdot_42, Sdot_52, Sdot_62;
      Scalar Sdot_13, Sdot_23, Sdot_33, Sdot_43, Sdot_53, Sdot_63;

      // Derivative of dndotXX_dqdot0 with respect to q0 and q1
      Scalar d_dndotx_dqdot1_dq1 = -s1; // dndotx_dqdot1 = c1;

      Scalar d_dndoty_dqdot0_dq0 = s0 * c1; // dndoty_dqdot0 = - c0 * c1;
      Scalar d_dndoty_dqdot0_dq1 = c0 * s1;

      // Scalar d_dndoty_dqdot1_dq0 = c0 * s1; // dndoty_dqdot1 = s0 * s1;
      Scalar d_dndoty_dqdot1_dq1 = s0 * c1;

      Scalar d_dndotz_dqdot0_dq0 = -c1 * c0; // dndotz_dqdot0 = - c1 * s0;
      Scalar d_dndotz_dqdot0_dq1 = s0 * s1;

      // Scalar d_dndotz_dqdot1_dq0 = s0 * s1; // dndotz_dqdot1 = - c0 * s1;
      Scalar d_dndotz_dqdot1_dq1 = -c0 * c1;

      // Upper part (translation)
      // Row 1, Column 1
      Sdot_11 =
        qdot0
          * (-dndoty_dqdot0 * radius_b * (-c0 * c2 * s1 + s0 * s2) + dndotz_dqdot0 * radius_c * (c0 * s2 + c2 * s0 * s1) + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq0 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq0)
        + qdot1
            * (dndoty_dqdot0 * radius_b * c1 * c2 * s0 - dndotz_dqdot0 * radius_c * c0 * c1 * c2 + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq1)
        - qdot2
            * (dndoty_dqdot0 * radius_b * (-c0 * c2 + s0 * s1 * s2) - dndotz_dqdot0 * radius_c * (c0 * s1 * s2 + c2 * s0));

      // Row 1, Column 2
      Sdot_12 =
        qdot0
          * (-dndoty_dqdot1 * radius_b * (-c0 * c2 * s1 + s0 * s2) + dndotz_dqdot1 * radius_c * (c0 * s2 + c2 * s0 * s1) + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq1)
        + qdot1 * (-dndotx_dqdot1 * radius_a * c2 * s1 + dndoty_dqdot1 * radius_b * c1 * c2 * s0 - dndotz_dqdot1 * radius_c * c0 * c1 * c2 + radius_a * c1 * c2 * d_dndotx_dqdot1_dq1 + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot1_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot1_dq1) - qdot2 * (dndotx_dqdot1 * radius_a * c1 * s2 + dndoty_dqdot1 * radius_b * (-c0 * c2 + s0 * s1 * s2) - dndotz_dqdot1 * radius_c * (c0 * s1 * s2 + c2 * s0));

      // Row 1, Column 3
      Sdot_13 = Scalar(0);

      // Row 2, Column 1
      Sdot_21 =
        -qdot0
          * (dndoty_dqdot0 * radius_b * (c0 * s1 * s2 + c2 * s0) + dndotz_dqdot0 * radius_c * (-c0 * c2 + s0 * s1 * s2) + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq0 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq0)
        - qdot1
            * (dndoty_dqdot0 * radius_b * c1 * s0 * s2 - dndotz_dqdot0 * radius_c * c0 * c1 * s2 + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq1 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq1)
        - qdot2
            * (dndoty_dqdot0 * radius_b * (c0 * s2 + c2 * s0 * s1) + dndotz_dqdot0 * radius_c * (-c0 * c2 * s1 + s0 * s2));

      // Row 2, Column 2
      Sdot_22 =
        -qdot0
          * (dndoty_dqdot1 * radius_b * (c0 * s1 * s2 + c2 * s0) + dndotz_dqdot1 * radius_c * (-c0 * c2 + s0 * s1 * s2) + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq1 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq1)
        + qdot1 * (dndotx_dqdot1 * radius_a * s1 * s2 - dndoty_dqdot1 * radius_b * c1 * s0 * s2 + dndotz_dqdot1 * radius_c * c0 * c1 * s2 - radius_a * c1 * s2 * d_dndotx_dqdot1_dq1 - radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot1_dq1 + radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot1_dq1) - qdot2 * (dndotx_dqdot1 * radius_a * c1 * c2 + dndoty_dqdot1 * radius_b * (c0 * s2 + c2 * s0 * s1) + dndotz_dqdot1 * radius_c * (-c0 * c2 * s1 + s0 * s2));

      // Row 2, Column 3
      Sdot_23 = Scalar(0);

      // Row 3, Column 1
      Sdot_31 =
        -qdot0 * c1
          * (dndoty_dqdot0 * radius_b * c0 + dndotz_dqdot0 * radius_c * s0 + radius_b * s0 * d_dndoty_dqdot0_dq0 - radius_c * c0 * d_dndotz_dqdot0_dq0)
        + qdot1
            * (-c1 * (radius_b * s0 * d_dndoty_dqdot0_dq1 - radius_c * c0 * d_dndotz_dqdot0_dq1) + s1 * (dndoty_dqdot0 * radius_b * s0 - dndotz_dqdot0 * radius_c * c0));

      // Row 3, Column 2
      Sdot_32 =
        -qdot0 * c1
          * (dndoty_dqdot1 * radius_b * c0 + dndotz_dqdot1 * radius_c * s0 + radius_b * s0 * d_dndoty_dqdot0_dq1 - radius_c * c0 * d_dndotz_dqdot0_dq1)
        + qdot1
            * (dndotx_dqdot1 * radius_a * c1 + dndoty_dqdot1 * radius_b * s0 * s1 - dndotz_dqdot1 * radius_c * c0 * s1 + radius_a * s1 * d_dndotx_dqdot1_dq1 - radius_b * c1 * s0 * d_dndoty_dqdot1_dq1 + radius_c * c0 * c1 * d_dndotz_dqdot1_dq1);

      // Row 3, Column 3
      Sdot_33 = Scalar(0);

      // Angular part (rows 4-6)
      Sdot_41 = -(qdot1 * c2 * s1 + qdot2 * c1 * s2);
      Sdot_51 = qdot1 * s1 * s2 - qdot2 * c1 * c2;
      Sdot_61 = qdot1 * c1;

      Sdot_42 = qdot2 * c2;
      Sdot_52 = -qdot2 * s2;
      Sdot_62 = Scalar(0);

      Sdot_43 = Scalar(0);
      Sdot_53 = Scalar(0);
      Sdot_63 = Scalar(0);

      data.Sdot.matrix() << Sdot_11, Sdot_12, Sdot_13, Sdot_21, Sdot_22, Sdot_23, Sdot_31, Sdot_32,
        Sdot_33, Sdot_41, Sdot_42, Sdot_43, Sdot_51, Sdot_52, Sdot_53, Sdot_61, Sdot_62, Sdot_63;
    }
    
    template<typename ConfigVector, typename TangentVector>
    void computeBiais(
      JointDataDerived & data,
      const Eigen::MatrixBase<ConfigVector> & qs,
      const Eigen::MatrixBase<TangentVector> & vs) const
    {
      Scalar c0, s0;
      SINCOS(qs(0), &s0, &c0);
      Scalar c1, s1;
      SINCOS(qs(1), &s1, &c1);
      Scalar c2, s2;
      SINCOS(qs(2), &s2, &c2);

      Scalar dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0, dndotz_dqdot1;
      dndotx_dqdot1 = c1;
      dndoty_dqdot0 = -c0 * c1;
      dndoty_dqdot1 = s0 * s1;
      dndotz_dqdot0 = -c1 * s0;
      dndotz_dqdot1 = -c0 * s1;

      computeBiais(
        s0, c0, s1, c1, s2, c2, dndotx_dqdot1, dndoty_dqdot0, dndoty_dqdot1, dndotz_dqdot0,
        dndotz_dqdot1, data);
    }

    void computeBiais(
      const Scalar & s0,
      const Scalar & c0,
      const Scalar & s1,
      const Scalar & c1,
      const Scalar & s2,
      const Scalar & c2,
      const Scalar & dndotx_dqdot1,
      const Scalar & dndoty_dqdot0,
      const Scalar & dndoty_dqdot1,
      const Scalar & dndotz_dqdot0,
      const Scalar & dndotz_dqdot1,
      JointDataDerived & data) const
    {
      // Compute Sdot for bias acceleration
      Scalar qdot0, qdot1, qdot2;
      qdot0 = data.joint_v(0);
      qdot1 = data.joint_v(1);
      qdot2 = data.joint_v(2);

      // last columns and last element of the second column are zero, 
      //  so we do not compute them
      Scalar Sdot_11, Sdot_21, Sdot_31, Sdot_41, Sdot_51, Sdot_61;
      Scalar Sdot_12, Sdot_22, Sdot_32, Sdot_42, Sdot_52;


      // Derivative of dndotXX_dqdot0 with respect to q0 and q1
      Scalar d_dndotx_dqdot1_dq1 = -s1; // dndotx_dqdot1 = c1;

      Scalar d_dndoty_dqdot0_dq0 = s0 * c1; // dndoty_dqdot0 = - c0 * c1;
      Scalar d_dndoty_dqdot0_dq1 = c0 * s1;

      // Scalar d_dndoty_dqdot1_dq0 = c0 * s1; // dndoty_dqdot1 = s0 * s1;
      Scalar d_dndoty_dqdot1_dq1 = s0 * c1;

      Scalar d_dndotz_dqdot0_dq0 = -c1 * c0; // dndotz_dqdot0 = - c1 * s0;
      Scalar d_dndotz_dqdot0_dq1 = s0 * s1;

      // Scalar d_dndotz_dqdot1_dq0 = s0 * s1; // dndotz_dqdot1 = - c0 * s1;
      Scalar d_dndotz_dqdot1_dq1 = -c0 * c1;

      // Upper part (translation)
      // Row 1, Column 1
      Sdot_11 =
        qdot0
          * (-dndoty_dqdot0 * radius_b * (-c0 * c2 * s1 + s0 * s2) + dndotz_dqdot0 * radius_c * (c0 * s2 + c2 * s0 * s1) + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq0 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq0)
        + qdot1
            * (dndoty_dqdot0 * radius_b * c1 * c2 * s0 - dndotz_dqdot0 * radius_c * c0 * c1 * c2 + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq1)
        - qdot2
            * (dndoty_dqdot0 * radius_b * (-c0 * c2 + s0 * s1 * s2) - dndotz_dqdot0 * radius_c * (c0 * s1 * s2 + c2 * s0));

      // Row 1, Column 2
      Sdot_12 =
        qdot0
          * (-dndoty_dqdot1 * radius_b * (-c0 * c2 * s1 + s0 * s2) + dndotz_dqdot1 * radius_c * (c0 * s2 + c2 * s0 * s1) + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot0_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot0_dq1)
        + qdot1 * (-dndotx_dqdot1 * radius_a * c2 * s1 + dndoty_dqdot1 * radius_b * c1 * c2 * s0 - dndotz_dqdot1 * radius_c * c0 * c1 * c2 + radius_a * c1 * c2 * d_dndotx_dqdot1_dq1 + radius_b * (c0 * s2 + c2 * s0 * s1) * d_dndoty_dqdot1_dq1 + radius_c * (-c0 * c2 * s1 + s0 * s2) * d_dndotz_dqdot1_dq1) - qdot2 * (dndotx_dqdot1 * radius_a * c1 * s2 + dndoty_dqdot1 * radius_b * (-c0 * c2 + s0 * s1 * s2) - dndotz_dqdot1 * radius_c * (c0 * s1 * s2 + c2 * s0));

      // Row 2, Column 1
      Sdot_21 =
        -qdot0
          * (dndoty_dqdot0 * radius_b * (c0 * s1 * s2 + c2 * s0) + dndotz_dqdot0 * radius_c * (-c0 * c2 + s0 * s1 * s2) + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq0 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq0)
        - qdot1
            * (dndoty_dqdot0 * radius_b * c1 * s0 * s2 - dndotz_dqdot0 * radius_c * c0 * c1 * s2 + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq1 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq1)
        - qdot2
            * (dndoty_dqdot0 * radius_b * (c0 * s2 + c2 * s0 * s1) + dndotz_dqdot0 * radius_c * (-c0 * c2 * s1 + s0 * s2));

      // Row 2, Column 2
      Sdot_22 =
        -qdot0
          * (dndoty_dqdot1 * radius_b * (c0 * s1 * s2 + c2 * s0) + dndotz_dqdot1 * radius_c * (-c0 * c2 + s0 * s1 * s2) + radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot0_dq1 - radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot0_dq1)
        + qdot1 * (dndotx_dqdot1 * radius_a * s1 * s2 - dndoty_dqdot1 * radius_b * c1 * s0 * s2 + dndotz_dqdot1 * radius_c * c0 * c1 * s2 - radius_a * c1 * s2 * d_dndotx_dqdot1_dq1 - radius_b * (-c0 * c2 + s0 * s1 * s2) * d_dndoty_dqdot1_dq1 + radius_c * (c0 * s1 * s2 + c2 * s0) * d_dndotz_dqdot1_dq1) - qdot2 * (dndotx_dqdot1 * radius_a * c1 * c2 + dndoty_dqdot1 * radius_b * (c0 * s2 + c2 * s0 * s1) + dndotz_dqdot1 * radius_c * (-c0 * c2 * s1 + s0 * s2));


      // Row 3, Column 1
      Sdot_31 =
        -qdot0 * c1
          * (dndoty_dqdot0 * radius_b * c0 + dndotz_dqdot0 * radius_c * s0 + radius_b * s0 * d_dndoty_dqdot0_dq0 - radius_c * c0 * d_dndotz_dqdot0_dq0)
        + qdot1
            * (-c1 * (radius_b * s0 * d_dndoty_dqdot0_dq1 - radius_c * c0 * d_dndotz_dqdot0_dq1) + s1 * (dndoty_dqdot0 * radius_b * s0 - dndotz_dqdot0 * radius_c * c0));

      // Row 3, Column 2
      Sdot_32 =
        -qdot0 * c1
          * (dndoty_dqdot1 * radius_b * c0 + dndotz_dqdot1 * radius_c * s0 + radius_b * s0 * d_dndoty_dqdot0_dq1 - radius_c * c0 * d_dndotz_dqdot0_dq1)
        + qdot1
            * (dndotx_dqdot1 * radius_a * c1 + dndoty_dqdot1 * radius_b * s0 * s1 - dndotz_dqdot1 * radius_c * c0 * s1 + radius_a * s1 * d_dndotx_dqdot1_dq1 - radius_b * c1 * s0 * d_dndoty_dqdot1_dq1 + radius_c * c0 * c1 * d_dndotz_dqdot1_dq1);

      // Angular part (rows 4-6)
      Sdot_41 = -(qdot1 * c2 * s1 + qdot2 * c1 * s2);
      Sdot_51 = qdot1 * s1 * s2 - qdot2 * c1 * c2;
      Sdot_61 = qdot1 * c1;

      Sdot_42 = qdot2 * c2;
      Sdot_52 = -qdot2 * s2;
      
      data.c.toVector() << Sdot_11 * qdot0 + Sdot_12 * qdot1,
                          Sdot_21 * qdot0 + Sdot_22 * qdot1,
                          Sdot_31 * qdot0 + Sdot_32 * qdot1,
                          Sdot_41 * qdot0 + Sdot_42 * qdot1,
                          Sdot_51 * qdot0 + Sdot_52 * qdot1,
                          Sdot_61 * qdot0;
    }
  }; // struct JointModelEllipsoidTpl

} // namespace pinocchio

#include <boost/type_traits.hpp>

namespace boost
{
  template<typename Scalar, int Options>
  struct has_nothrow_constructor<::pinocchio::JointModelEllipsoidTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_copy<::pinocchio::JointModelEllipsoidTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_constructor<::pinocchio::JointDataEllipsoidTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };

  template<typename Scalar, int Options>
  struct has_nothrow_copy<::pinocchio::JointDataEllipsoidTpl<Scalar, Options>>
  : public integral_constant<bool, true>
  {
  };
} // namespace boost

#endif // ifndef __pinocchio_multibody_joint_spherical_ZYX_hpp__
