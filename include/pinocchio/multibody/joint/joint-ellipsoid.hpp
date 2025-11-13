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

  /// \brief Ellipsoid joint - constrains motion to ellipsoid surface with 3-DOF.
  ///
  /// The configuration space uses three angles (q₀, q₁, q₂) representing:
  /// - Rotation about the x-axis
  /// - Rotation about the y-axis
  /// - Spin about the "normal" direction
  ///
  /// The joint position on the ellipsoid surface is computed as:
  /// \f$ \mathbf{p} = (a \sin q_1, -b \sin q_0 \cos q_1, c \cos q_0 \cos q_1) \f$
  ///
  /// where \f$ a, b, c \f$ are the radii along the x, y, z axes respectively.
  ///
  /// \note For non-spherical ellipsoids, the third rotation axis is only approximately
  /// normal to the surface. It corresponds to the normal of an equivalent sphere while
  /// the translation follows the true ellipsoid surface. The "normal" direction is
  /// truly normal only when all radii are equal (sphere case).
  ///
  /// \sa Seth et al., "Minimal formulation of joint motion for biomechanisms,"
  /// Nonlinear Dynamics 62(1):291-303, 2010.
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

    /// @brief Default constructor. Creates a sphere (0.01, 0.01, 0.01).
    JointModelEllipsoidTpl()
    : radius_a(Scalar(0.01))
    , radius_b(Scalar(0.01))
    , radius_c(Scalar(0.01))
    {
    }

    /// @brief Constructor with specified radii.
    /// @param a Semi-axis length along x-direction
    /// @param b Semi-axis length along y-direction
    /// @param c Semi-axis length along z-direction
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

    /// @brief Computes the spatial transformation M(q) from joint configuration.
    /// @param[in] c0, s0 Cosine and sine of q[0]
    /// @param[in] c1, s1 Cosine and sine of q[1]
    /// @param[in] c2, s2 Cosine and sine of q[2]
    /// @param[out] data Joint data where M will be stored
    void computeSpatialTransform(
      Scalar c0, Scalar s0, Scalar c1, Scalar s1, Scalar c2, Scalar s2, JointDataDerived & data)
      const
    {
      // clang-format off
      data.M.rotation() << c1 * c2                , -c1 * s2                , s1      ,
                           c0 * s2 + c2 * s0 * s1 , c0 * c2 - s0 * s1 * s2  ,-c1 * s0 ,
                          -c0 * c2 * s1 + s0 * s2 , c0 * s1 * s2 + c2 * s0  , c0 * c1;
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

      computeSpatialTransform(c0, s0, c1, s1, c2, s2, data);

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

      data.S.matrix() << S_11   , S_12  , Scalar(0),
                         S_21   , S_22  , Scalar(0),
                         S_31   , S_32  , Scalar(0),
                         c1 * c2, s2    , Scalar(0),
                        -c1 * s2, c2    , Scalar(0),
                         s1     , Scalar(0), Scalar(1);
      // clang-format on
    }

    /// @brief Computes the bias acceleration c(q, v) = Sdot(q)·v.
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
      // so we do not compute them
      Scalar Sdot_11, Sdot_21, Sdot_31, Sdot_41, Sdot_51, Sdot_61;
      Scalar Sdot_12, Sdot_22, Sdot_32, Sdot_42, Sdot_52;

      // Derivative of dndotXX_dqdot0 with respect to q0 and q1
      Scalar d_dndotx_dqdot1_dq1 = -s1;

      Scalar d_dndoty_dqdot0_dq0 = s0 * c1;
      Scalar d_dndoty_dqdot0_dq1 = c0 * s1;

      Scalar d_dndoty_dqdot1_dq1 = s0 * c1;

      Scalar d_dndotz_dqdot0_dq0 = -c1 * c0;
      Scalar d_dndotz_dqdot0_dq1 = s0 * s1;

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

      data.c.toVector() << Sdot_11 * qdot0 + Sdot_12 * qdot1, Sdot_21 * qdot0 + Sdot_22 * qdot1,
        Sdot_31 * qdot0 + Sdot_32 * qdot1, Sdot_41 * qdot0 + Sdot_42 * qdot1,
        Sdot_51 * qdot0 + Sdot_52 * qdot1, Sdot_61 * qdot0;
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

#endif // ifndef __pinocchio_multibody_joint_ellipsoid_hpp__
