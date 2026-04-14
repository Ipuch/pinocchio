//
// Copyright (c) 2024-2025 INRIA
//

#pragma once

// IWYU pragma: private, include "pinocchio/constraints.hpp"

#ifdef PINOCCHIO_LSP
  #undef PINOCCHIO_LSP
  #include "pinocchio/constraints.hpp"
#endif // PINOCCHIO_LSP

namespace pinocchio
{
  template<typename _Scalar, int _Options>
  struct ConstraintCollectionDefaultTpl
  {
    typedef _Scalar Scalar;
    static constexpr int Options = _Options;

    typedef PointAnchorConstraintModelTpl<Scalar, Options> PointAnchorConstraintModel;
    typedef PointAnchorConstraintDataTpl<Scalar, Options> PointAnchorConstraintData;

    typedef PointContactConstraintModelTpl<Scalar, Options> PointContactConstraintModel;
    typedef PointContactConstraintDataTpl<Scalar, Options> PointContactConstraintData;

    typedef JointFrictionConstraintModelTpl<Scalar, Options> JointFrictionConstraintModel;
    typedef JointFrictionConstraintDataTpl<Scalar, Options> JointFrictionConstraintData;

    typedef JointLimitConstraintModelTpl<Scalar, Options> JointLimitConstraintModel;
    typedef JointLimitConstraintDataTpl<Scalar, Options> JointLimitConstraintData;

    typedef CoordinateCoupletConstraintModelTpl<Scalar, Options> CoordinateCoupletConstraintModel;
    typedef CoordinateCoupletConstraintDataTpl<Scalar, Options> CoordinateCoupletConstraintData;

    typedef FrameAnchorConstraintModelTpl<Scalar, Options> FrameAnchorConstraintModel;
    typedef FrameAnchorConstraintDataTpl<Scalar, Options> FrameAnchorConstraintData;

    typedef boost::variant<
      BlankConstraintModel,
      PointAnchorConstraintModel,
      PointContactConstraintModel,
      JointFrictionConstraintModel,
      JointLimitConstraintModel,
      CoordinateCoupletConstraintModel,
      FrameAnchorConstraintModel>
      ConstraintModelVariant;

    typedef boost::variant<
      BlankConstraintData,
      PointAnchorConstraintData,
      PointContactConstraintData,
      JointFrictionConstraintData,
      JointLimitConstraintData,
      CoordinateCoupletConstraintData,
      FrameAnchorConstraintData>
      ConstraintDataVariant;
  }; // struct ConstraintCollectionDefaultTpl

  typedef ConstraintCollectionDefault::ConstraintModelVariant ConstraintModelVariant;
  typedef ConstraintCollectionDefault::ConstraintDataVariant ConstraintDataVariant;

} // namespace pinocchio
