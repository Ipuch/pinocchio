//
// Copyright (c) 2026 INRIA
// Copyright (c) 2026 ISIR
//
#include "pinocchio/bindings/python/algorithm/algorithms.hpp"
#include "pinocchio/algorithm/spline-regression.hpp"

namespace pinocchio
{
  namespace python
  {

    static context::JointModelSpline fitSplineJoint_proxy(
      const std::vector<context::SE3> & data_frames,
      const context::VectorXs & q_data,
      const size_t nb_control_frames,
      const size_t degree,
      SplineRegressionSettings & settings)
    {
      return fitSplineJoint(data_frames, q_data, nb_control_frames, degree, settings);
    }

    static context::JointModelSpline fitSplineJoint_default_proxy(
      const std::vector<context::SE3> & data_frames,
      const context::VectorXs & q_data,
      const size_t nb_control_frames,
      const size_t degree)
    {
      return fitSplineJoint(data_frames, q_data, nb_control_frames, degree);
    }

    void exposeSplineRegression()
    {
      typedef context::Scalar Scalar;

      bp::class_<SplineRegressionSettings>(
        "SplineRegressionSettings",
        "Settings for the spline joint regression algorithm (see fitSplineJoint).",
        bp::init<>(bp::args("self"), "Default constructor."))
        .def(bp::init<Scalar, Scalar, int>(
          bp::args("self", "accuracy", "damping", "max_iter"),
          "Constructor with the main setting parameters."))
        .def_readwrite(
          "absolute_accuracy", &SplineRegressionSettings::absolute_accuracy,
          "Absolute accuracy: stop when the residual goes below this value.")
        .def_readwrite(
          "relative_accuracy", &SplineRegressionSettings::relative_accuracy,
          "Relative accuracy: stop when the relative residual improvement between two "
          "iterates goes below this value.")
        .def_readwrite(
          "damping", &SplineRegressionSettings::damping,
          "Tikhonov damping added to the normal matrix (useful on noisy or sparse data).")
        .def_readwrite(
          "smoothing_weight", &SplineRegressionSettings::smoothing_weight,
          "Hessian smoothing weight in [0, 1) (0 disables smoothing).")
        .def_readwrite(
          "smoothing_drag", &SplineRegressionSettings::smoothing_drag,
          "Hessian smoothing drag rate in (0, 1].")
        .def_readwrite(
          "max_iter", &SplineRegressionSettings::max_iter,
          "Maximal number of Gauss-Newton iterations.")
        .def_readwrite(
          "absolute_residual", &SplineRegressionSettings::absolute_residual,
          "Final data residual (output).")
        .def_readwrite(
          "relative_residual", &SplineRegressionSettings::relative_residual,
          "Relative residual improvement at the last iterate (output).")
        .def_readwrite(
          "iter", &SplineRegressionSettings::iter,
          "Total number of iterations performed (output).");

      bp::def(
        "fitSplineJoint", fitSplineJoint_proxy,
        bp::args("data_frames", "q_data", "nb_control_frames", "degree", "settings"),
        "Fit the control frames of a spline joint to data placements by Gauss-Newton "
        "regression in se(3) (Lee & Terzopoulos 2008, Sec. 6.2-6.3).\n\n"
        "The returned joint uses an open uniform knot vector spanning the q_data range, "
        "so it is directly driven by the data parameter (e.g. a physical joint angle).\n\n"
        "Parameters:\n"
        "\tdata_frames: list of SE3 data placements to fit (size N)\n"
        "\tq_data: sample parameters, strictly increasing (size N)\n"
        "\tnb_control_frames: number of spline control frames (> degree, <= N)\n"
        "\tdegree: degree of the B-spline basis functions\n"
        "\tsettings: regression settings; on output it contains the final residual and "
        "the number of iterations\n");

      bp::def(
        "fitSplineJoint", fitSplineJoint_default_proxy,
        bp::args("data_frames", "q_data", "nb_control_frames", "degree"),
        "Fit the control frames of a spline joint to data placements by Gauss-Newton "
        "regression in se(3), with default regression settings.\n");
    }

  } // namespace python
} // namespace pinocchio
