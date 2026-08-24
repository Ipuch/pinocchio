---
name: pinocchio-constraint
description: Add a new constraint model + data to Pinocchio's C++ constraint API (the ConstraintModelTpl/ConstraintDataTpl variants) — distance, surface, coupling, joint-space constraints, etc. Covers which CRTP base to derive from, the full registration checklist (a missed file is a hard compile error), the Eigen/traits traps for residual sizes other than 3 and 6, the invariants worth testing, and a fast single-TU compile-and-run loop. Use when asked to implement/add/extend a pinocchio constraint or "contrainte".
---

# Building a constraint in Pinocchio

## Read first

- `development/convention.md` — non-negotiable. Public `.hpp` in `include/pinocchio/`, private
  `.hxx` in `include/pinocchio/src/`. **A private header may not include anything**; the omnibus
  public header (`constraints.hpp`) pulls in every dependency in the right order. Every `.hxx`
  starts with `#pragma once`, the `// IWYU pragma: private, include ...` line, and the
  `PINOCCHIO_LSP` guard block.
- The nearest existing constraint. `point-anchor-constraint.hxx` (3D, trivial), 
  `point-contact-constraint.hxx` (3D + a stored parameter), `joint-limit-constraint.hxx`
  (dynamic size, hardest), `constant-length-constraint.hxx` (scalar — **the reference for any
  size other than 3 or 6**).

All paths below are relative to `include/pinocchio/src/constraints/` unless stated otherwise.

## Step 1 — pick the CRTP base

| Constraint depends on | Residual | Derive from | You must write |
|---|---|---|---|
| relative position of two points | 3 | `PointConstraintModelBase` | traits + ctors + `setImpl`; the algebra is inherited |
| relative placement of two frames | 6 | `FrameConstraintModelBase` | same |
| joint coordinates (limits, friction) | dynamic | `JointWiseConstraintModelBase` | a lot; read `joint-limit-constraint.hxx` |
| **anything else** (scalar, 2, 4, …) | any | `BinaryKinematicsConstraintModelBase` | the **whole** traits block + ~12 `*Impl` methods |

The last row is the interesting one, and it is cheaper than it looks whenever the new residual is a
smooth scalar/low-rank function of a *point* constraint:

> If `φ(q) = f(x(q))` with `x` the relative point position, then `J_φ = ∇f(x)ᵀ · ∂x/∂q`, and the
> whole Pinocchio projector machinery follows: `A1 = ∇f(x)ᵀ · A1_point`, `A2 = ∇f(x)ᵀ · A2_point`.

The 3×6 projectors of the underlying point constraint already exist as
`internal::relativePointProjector1/2` in `relative-point-projectors.hxx` — call them and
left-multiply by the gradient row, do not re-derive them. Two constraints already do this:

| constraint | `f(x)` | `∇f` |
|---|---|---|
| `constant-length-constraint.hxx` | `‖x‖ − L` | `u = x/‖x‖` |
| `ellipsoid-point-constraint.hxx` | `(xᵀAx − 1) / (2‖Ax‖)`, `A = diag(1/aᵢ²)` | `Ax/‖Ax‖ − φ_alg·A²x/(2‖Ax‖³)` |

Copy whichever is closer, swap `f`, and you are most of the way there. The acceleration error
needs the **Hessian** too — write `ẋᵀ∇²f ẋ` as scalar contractions of `ẋ` with `Ax`, `A²x`, … and
let test 4 below (FD + `J·a + err(a=0) == err(a)`) tell you whether you got it right.

Two things worth stating explicitly, because they are what makes such a residual usable:

- **Scale the residual to metres.** For a level set `g(x) = 0`, use `g/‖∇g‖` rather than `g`: it is
  the signed distance to first order, so it is comparable with `eps`, with marker errors, and with
  the other terms of a merit function. On the manifold `‖∇φ‖ == 1` exactly — a free assertion.
- **Differentiate the scaling factor too.** Freezing `1/(2‖Ax‖)` inside `calc` gives a gradient
  that is right *on* the manifold and a few percent off a few millimetres away — fine for a Newton
  step, wrong for a line search that starts far from it. Test the Jacobian against finite
  differences at several distances from the surface, not only on it.

## Step 2 — the header

Contract the base classes expect (`constraint-model-base.hxx` lists it at the bottom):
`classnameImpl`, `shortnameImpl`, `createDataImpl`, `setImpl`, `calcImpl`, `jacobianImpl`,
`jacobianMatrixProductImpl` ×2, `jacobianTransposeMatrixProductImpl` ×2,
`mapConstraintForceToJointForcesImpl`, `mapJointMotionsToConstraintMotionImpl`,
`appendCouplingConstraintInertiasImpl` ×2 (vector + `MatrixBlockElementTpl`), `getA1Impl`,
`getA2Impl`, plus `cast`, `operator==/!=`.

`calcImpl` must fill, in cdata: `oMc1`, `oMc2`, `c1Mc2`, `constraint_position_error`,
`constraint_velocity_error`, `constraint_acceleration_error`, and the six `A{1,2,}_{world,local}`.
The position/velocity/acceleration triple must be consistent: velocity error `= J·v`, acceleration
error affine in `a` with `J` as its linear part. For `φ = f(x)`:

```
φ̇  = ∇f·ẋ
φ̈  = ∇f·ẍ + ẋᵀ ∇²f ẋ          (for f = ‖x‖ : φ̈ = u·ẍ + (‖ẋ‖² − (u·ẋ)²)/‖x‖, u = x/‖x‖)
```

Copy `ẋ`/`ẍ` verbatim from `PointConstraintModelBase::calcImpl` — the Coriolis terms are easy to
get wrong. Subtract `desired_constraint_offset/velocity/acceleration` like the point base does.

Store a **dedicated parameter** as a protected `m_foo` with `getFoo()/setFoo()` accessors
(like `PointContactConstraintModel::m_friction`) rather than overloading
`desired_constraint_offset`, and validate it in the setter with
`PINOCCHIO_CHECK_INPUT_ARGUMENT(check_expression_if_real<Scalar>(...), "msg")`.

Runtime guards inside `calc` are `assert(check_expression_if_real<Scalar>(...))`, never a throw:
`calc` is a hot loop, and `check_expression_if_real` is what keeps CasADi/CppAD scalars compiling.

`check_expression_if_real` is not enough on its own: **the expression inside it must still be a
scalar comparison**. `radii.minCoeff() > Scalar(0)` compiles for `double` and fails for
`casadi::SX` with "cannot convert `casadi::Matrix<casadi::SXElem>` to `bool`", because the Eigen
reduction needs `operator<` to return a bool. Loop over the components instead — which is why no
other constraint header contains a `minCoeff`/`maxCoeff`. Nothing in a default build catches this:
you need `-DBUILD_WITH_CASADI_SUPPORT=ON`, or the standalone probe below.

## Step 3 — registration checklist

Every line is a hard compile error (or a silently missing binding) if skipped.

| File | What to add |
|---|---|
| `src/constraints/<name>-constraint.hxx` | the new model + data |
| `src/constraints/fwd.hxx` | fwd decls + the two `context::Scalar` typedefs |
| `constraints.hpp` | `#include` the new `.hxx`, **after** the base it derives from |
| `src/constraints/constraint-collection-default.hxx` | member typedefs **and** both `boost::variant` lists |
| `src/constraints/utils.hxx` | `ComputeBlockDiagonalPatternImpl` specialization |
| `src/serialization/constraints-model.hxx` | `serialize()`; reach private members via an `internal::…Accessor` struct |
| `src/serialization/constraints-data.hxx` | `serialize()` for the data |
| `bindings/python/context/generic.hpp` | the 4 `context::` typedefs (model/data + vectors) |
| `bindings/python/algorithm/constraints/constraints-models.hpp` | `expose_constraint_model<>` — only if extra ctors/accessors |
| `bindings/python/algorithm/constraints/constraints-datas.hpp` | `expose_constraint_data<>` — **required** unless the data derives from `PointConstraintDataBase`/`FrameConstraintDataBase` |
| `unittest/serialization-constraints.cpp` | `initConstraint<>` specialization |
| `unittest/CMakeLists.txt` | `add_pinocchio_unit_test(<name>-constraint)` |
| `CHANGELOG.md` | entry under `## [Unreleased]` / `### Added` |

Python class registration itself is automatic: `expose-constraints.cpp` runs
`boost::mpl::for_each` over the variant types.

## Pitfalls (all verified the hard way)

- **Eigen rejects `Matrix<Scalar, 1, N, ColMajor>`** (`INVALID_MATRIX_TEMPLATE_PARAMETERS`) for
  `N != 1`. A `Size == 1` constraint must spell `Eigen::RowMajor` in `JacobianMatrixType`, in
  `JacobianMatrixProductReturnType`, and in the data's `A*` members. `Matrix<Scalar,1,1,Options>`
  is fine. `BinaryKinematicsConstraintModelBase::MatrixSize6` already switches storage order on
  `Size == 1` — keep it that way.
- **Traits are not inherited** when you derive straight from `BinaryKinematicsConstraintModelBase`
  (there is no `traits<BinaryKinematicsConstraintModelBase<D>>`). Your `traits<>` must define
  everything the point/frame traits do, *including* the nested `JacobianMatrixProductReturnType`
  and `JacobianTransposeMatrixProductReturnType` templates.
- `ComputeBlockDiagonalPatternImpl`'s primary template is a `static_assert(false)` — but it is
  `#if`-ed out on GCC < 13, where a missing specialization degrades to an empty block list instead
  of an error. Do not rely on the compiler to remind you. A bilateral equality constraint emits
  `{MatrixBlockType::ScalarIdentity, Size}`.
- `unittest/serialization-constraints.cpp` iterates the variant, so adding a type there breaks the
  build until `initConstraint<YourModel>` exists.
- In tests, **never name an `SE3` `M_E`**: `<cmath>` defines `M_E` as Euler's number and the error
  messages point at your constructor, not at the macro.
- `SE3::act()` wants a materialised `Vector3`, not an Eigen expression — `M.act(a.cwiseProduct(b))`
  fails with "has no member named `se3Action`". Assign to a `Vector3d` first.
- In `expose_constraint_model`, a setter taking a template argument needs the raw `+[](...)`
  lambda-to-function-pointer form when it is followed by `bp::args(...)`;
  `bp::make_function(+[]...)` only works when the docstring comes straight after.
- The Python **data** inheritance visitor is SFINAE-gated on `PointConstraintDataBase` /
  `FrameConstraintDataBase`. A new data family silently gets zero properties — add the
  `expose_constraint_data` specialization.
- `A1_world + A2_world == 0` (exactly) for any constraint invariant under a rigid displacement of
  the whole system. Cheap, sharp self-check — and it is why the "dof drives both joints" branch of
  `jacobianMatrixProduct` correctly contributes nothing.

## Step 4 — tests worth writing

Model the file on `unittest/point-anchor-constraint.cpp`. The checks that actually caught mistakes:

1. residual vanishes on the manifold; sign convention on either side of it;
2. `getA1/getA2` equal the gradient row times the equivalent `PointAnchorConstraintModel`
   projectors, in **both** `WorldFrameTag` and `LocalFrameTag` — this is the strongest single test;
3. `J == A1·J1 + A2·J2` in both frames, and `J` vs finite differences on `pin::integrate`;
4. `constraint_velocity_error == J·v`; velocity and acceleration errors vs finite differences;
   `J·a + error(a=0) == error(a)`;
5. `check_jacobians_operations(...)` from `unittest/constraints/jacobians-checker.hpp` — covers
   the four product overloads and `SetTo/AddTo/RmTo`;
6. Jacobian sparsity vs `getRowIndexes`; invariance under a rigid displacement of the base;
7. `cast`, compliance round-trip, the `ConstraintModel` variant path, and a Cholesky assembly
   against a hand-built KKT matrix.

## Step 5 — verify

```bash
# single translation unit, asserts on; ~40 s instead of a full rebuild
.claude/skills/pinocchio-constraint/scripts/build-test.sh unittest/<name>-constraint.cpp
.claude/skills/pinocchio-constraint/scripts/build-test.sh unittest/serialization-constraints.cpp
```

Then, because adding a variant type recompiles a lot of the library, syntax-check the heavy
consumers (run them in parallel, each takes a minute or two):

```
unittest/{constraint-variants,constraint-cholesky,loop-constrained-aba,delassus-operator-rigid-body}.cpp
unittest/{contact-dynamics,contact-inverse-dynamics,closed-loop-dynamics,constraint-jacobian}.cpp
src/algorithm/constraints/utils.cpp   src/algorithm/constraint-cholesky.cpp
```

The project builds with `-Wconversion -Wcast-qual -Wcast-align -Wwrite-strings -pedantic`; the
script inherits them from `compile_commands.json`, so **aim for zero warnings**, and re-run once
with `--release` since asserts and `-O3` change what gets instantiated.

**Check the symbolic scalars too**, without waiting for a CasADi-enabled build: write a tiny
`main()` that instantiates your model+data with `casadi::SX` and calls `calc` and `jacobian` on
them, then `-fsyntax-only` it against a prefix that has casadi (`.pixi/envs/all` or `.pixi/envs/casadi`).
Template instantiation is what surfaces the bool-conversion failures, so the probe must actually
call the methods, not merely name the type.

```cpp
typedef casadi::SX ADScalar;
pinocchio::ModelTpl<ADScalar> ad_model = model.cast<ADScalar>();
typename pinocchio::ModelTpl<ADScalar>::Data ad_data(ad_model);
// ... forwardKinematics + computeJointJacobians on a symbolic q ...
CM cmodel(ad_model, /* ... */);  typename CM::ConstraintData cdata(cmodel);
cmodel.calc(ad_model, ad_data, cdata);
const auto J = cmodel.jacobian(ad_model, ad_data, cdata);
```

Python bindings are often `OFF` in a local build tree. Type-check them directly:

```bash
NPINC=$(.pixi/envs/default/bin/python -c 'import numpy; print(numpy.get_include())')
.pixi/envs/default/bin/x86_64-conda-linux-gnu-c++ -std=c++17 \
  -DBOOST_MPL_LIMIT_LIST_SIZE=30 -DBOOST_MPL_LIMIT_VECTOR_SIZE=30 \
  -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_FUSION_INVOKE_MAX_ARITY=12 \
  -I build -I build/include -I include -I include/pinocchio/deprecated \
  -isystem .pixi/envs/default/include/eigen3 -isystem .pixi/envs/default/include \
  -isystem .pixi/envs/default/include/python3.* -isystem "$NPINC" \
  -fsyntax-only bindings/python/algorithm/constraints/expose-constraints.cpp
```

Finish with formatting — the PR template requires it:

```bash
pixi run lint                       # = pre-commit run --all
# offline fallback: use the clang-format pinned in .pre-commit-config.yaml, e.g.
#   find ~/.cache/pre-commit -name clang-format -type f   # then check --version matches the rev
```

## Step 6 — exercise it from Python, in a dedicated conda env

The unit tests prove the maths. This step proves the constraint is *usable*: a full build with
Python bindings, a worked example, and a Viser visualization.

### 6a. The environment

```bash
conda create -y -n pin4-dev-constraint -c conda-forge \
  python=3.12 cxx-compiler cmake ninja pkg-config ccache \
  libboost-devel libboost-python-devel eigen numpy scipy eigenpy urdfdom tinyxml tinyxml2 \
  coal meshcat-python matplotlib viser example-robot-data
```

> **Trap, hit every time on this machine.** A `pin_dev` entry is hard-wired at the front of `PATH`
> by the shell profile, and neither `conda activate` nor `conda run` displaces it — a configure
> that *looked* fine silently picked `cmake` and `python` from the wrong env. Always pass explicit
> paths, and **read back the configure log** to confirm Python/eigenpy/coal resolve inside
> `$P` before starting a 20-minute build.

### 6b. Configure and build

```bash
P=$HOME/miniconda3/envs/pin4-dev-constraint
export PATH=$P/bin:$PATH CMAKE_PREFIX_PATH=$P \
       CC=$P/bin/x86_64-conda-linux-gnu-cc CXX=$P/bin/x86_64-conda-linux-gnu-c++

$P/bin/cmake -G Ninja -B build_constraint -S . \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$P \
  -DPython3_EXECUTABLE=$P/bin/python -DPYTHON_EXECUTABLE=$P/bin/python \
  -DBUILD_PYTHON_INTERFACE=ON -DBUILD_WITH_COLLISION_SUPPORT=ON -DBUILD_WITH_URDF_SUPPORT=ON \
  -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARK=OFF -DBUILD_UTILS=OFF \
  -DBUILD_WITH_CASADI_SUPPORT=OFF -DBUILD_WITH_AUTODIFF_SUPPORT=OFF \
  -DBUILD_WITH_CODEGEN_SUPPORT=OFF -DBUILD_WITH_SDF_SUPPORT=OFF -DBUILD_WITH_OPENMP_SUPPORT=OFF \
  -DGENERATE_PYTHON_STUBS=OFF -DINSTALL_DOCUMENTATION=OFF \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

$P/bin/cmake --build build_constraint --parallel 20
```

- **Collision support is mandatory** to display anything: the visualizers take `coal` shapes.
- Use a build directory of its own. `build/` is the C++-only tree that Step 5's script targets;
  do not repoint it.
- `-DENABLE_TEMPLATE_INSTANTIATION=OFF` speeds up building *one* C++ test, but slows the bindings
  down — leave it ON here.

**Install into the env** — otherwise `conda activate pin4-dev-constraint; python -c "import
pinocchio"` fails with `ModuleNotFoundError`, because building alone copies nothing into
`$CONDA_PREFIX`:

```bash
$P/bin/cmake --install build_constraint      # `make install` in the recipe
```

While iterating you can skip the install and point at the build tree instead — but then say so,
and remember the env itself is still empty:

```bash
PYTHONPATH=$PWD/build_constraint/bindings/python $P/bin/python examples/<name>.py
```

### 6c. The example

Put it in `examples/<constraint-name>-kinematics.py` and register it in
`examples/CMakeLists.txt` — in the `BUILD_WITH_COLLISION_SUPPORT` block if it does `import coal`
at module level, in the plain `${PROJECT_NAME}_PYTHON_EXAMPLES` list otherwise. Examples run as
CI tests, so it must terminate and must never block.

**Match the house style**: these examples are *linear top-level scripts*. No `main()`, no
`argparse`, no helper functions — constants, then the model built statement by statement, then
the computation, then the viewer. Read `examples/simulation-closed-kinematic-chains.py` first and
copy its shape, including `pin.XAxis`/`pin.ZAxis`, `pin.Quaternion.FromTwoVectors`,
`pin.Inertia.FromBox`, and `import coal` (not `pin.coal`).

The natural demo for a kinematic constraint is **constrained forward kinematics**: drive some
dofs, solve for the rest with a Newton iteration on the residual.

```python
pin.forwardKinematics(model, data, q, zero, zero)   # calc() reads data.v and data.a
pin.computeJointJacobians(model, data, q)           # jacobian() needs this on top
cmodel.calc(model, data, cdata)
residual = cdata.constraint_position_error[0]
J = cmodel.jacobian(model, data, cdata)[:, free_idx]
dq_free = np.linalg.lstsq(J, -np.array([residual]), rcond=None)[0]
q = pin.integrate(model, q, dq)
```

Assert the final residual (`< 1e-10 m`) so the example is a real regression test, and show that
`J @ v == 0` characterises the admissible velocities.

For a **dynamics** example, note that `pin.constraintDynamics` and `pin.lcaba` are exposed in
Python for `RigidConstraintModel` only — the new constraint types are commented out in
`bindings/python/algorithm/expose-constrained-dynamics.cpp`. Assemble the KKT system yourself; it
is a handful of lines and it exercises all three quantities `calc()` produces:

```python
pin.forwardKinematics(model, data, q, v, np.zeros(model.nv))   # a = 0 -> gamma
pin.computeJointJacobians(model, data, q)
cmodel.calc(model, data, cdata)
phi   = cdata.constraint_position_error[0]
phid  = cdata.constraint_velocity_error[0]
gamma = cdata.constraint_acceleration_error[0]     # drift, since a = 0
# [ M  -J^T ] [ddq   ]   [tau - b                   ]
# [ J    0  ] [lambda] = [-gamma - Kd*phid - Kp*phi ]
```

Validate it against something analytic — a spherical pendulum's tension
`m (g cos θ + v² / L)` matched to 4 digits is a far stronger check than "it looks right".
Baumgarte gains live on the model (`cmodel.baumgarte_corrector_parameters.Kp/.Kd`) but `calc`
does not apply them: the algorithm does. With `Kp = 1e4`, `Kd = 2 sqrt(Kp)` and `dt = 1e-3`,
expect a residual around `1e-4 m` — set the assertion from a measured run, not from hope.

Two binding traps, both found by running the example rather than by reading the code:

- **A scalar constraint's `jacobian()` comes back 1-D.** `JacobianMatrixType` is a `1 x nv` Eigen
  *row vector*, and eigenpy maps every Eigen vector to a flat numpy array — so you get shape
  `(nv,)`, not `(1, nv)` as for a 3D constraint. Wrap it in `np.atleast_2d(...)` before slicing
  columns. Same story for `constraint_position_error`: a 1-element array, hence the `[0]`.
- **`pin.GeometryObject(name, parent_joint, placement, geometry)`** — placement *before* geometry.
  Some examples in `examples/` still use the opposite order and raise `Boost.Python.ArgumentError`.

### 6d. Viser

`from pinocchio.visualize import ViserVisualizer`; `viz.viewer` is the raw `viser.ViserServer`,
so anything the model does not own — and a constraint owns nothing — is drawn through it:

```python
viz.viewer.scene.add_line_segments(
    "/coupler",
    points=np.stack([p1, p2])[None].astype(np.float32),  # (n_segments, 2, 3)
    colors=(240, 190, 40),
    thickness=0.02,                                      # NOT line_width
)
```

Verified against viser 1.1: handles are **immutable** (`points` has no setter), so animate by
re-adding under the same name each frame — viser replaces the node.

Two more things that only show up at runtime:

- **`initViewer(open=True)` blocks forever** on `while len(viewer.get_clients()) == 0` when no
  browser connects. Fine when you run it by hand, fatal in CI. Use `open=False` and print
  `http://{viz.viewer.get_host()}:{viz.viewer.get_port()}` instead. (Meshcat's `open=True` does
  not block, which is why the Meshcat examples get away with it.)
- **`ViserVisualizer` handles Box, Sphere, Cylinder, Convex and meshes — nothing else**, and
  raises `RuntimeError: Unsupported geometry type` for the rest. Capsule and Ellipsoid both fall
  in that hole, although the Meshcat examples use capsules freely. Swap a capsule for a cylinder;
  for an ellipsoid, scale a `trimesh.creation.icosphere` by the radii and push it through
  `viz.viewer.scene.add_mesh_simple(name, vertices, faces, color=..., opacity=...)` — which is the
  right thing anyway, since a constraint manifold is not a body of the model.

Check the signature before writing the call — this API is not stable across viser versions:

```bash
$P/bin/python -c "import viser,inspect; print(inspect.signature(viser.SceneApi.add_line_segments))"
```
