#!/usr/bin/env python3
"""DOLFINx reference for associated, perfectly plastic Drucker--Prager.

This is a plane-strain strip-footing benchmark.  The constitutive equations and
benchmark parameters follow Cermak, Sysala, and Valdman (2019):

  https://doi.org/10.1016/j.amc.2019.02.054
  https://github.com/matlabfem/matlab_fem_elastoplasticity

Sign convention: tension is positive.  Symmetric tensors use the Mandel vector

    [xx, yy, zz, sqrt(2) xy],

so tensor inner products are ordinary vector dot products.  The yield function
is

    f(sigma) = sqrt(J2) + eta p - c <= 0,    p = tr(sigma) / 3.

The friction-angle and cohesion conversion is the plane-strain approximation
used by the reference implementation.  The flow rule is associated and there
is no hardening or tension cut-off.

Run the inexpensive constitutive checks only with

    python3 drucker_prager_fenicsx.py --self-test-only

or run the DOLFINx strip-footing problem and save its pressure/settlement path:

    python3 drucker_prager_fenicsx.py --output drucker_prager_results.csv
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# ==============================================================================
#  ------------------------------------------------------------
#  ------- Algorithm below is explained here -------
#  ------------------------------------------------------------
#    run_constitutive_self_tests()
#      check_elastic_response()
#      check_smooth_cone_return()
#      check_apex_return()
#      check_consistent_tangent_by_finite_difference()
#    Loop 1 : prescribed settlement increments       # Load loop
#      Loop 2 : max_newton_iterations                 # Newton loop
#        evaluate_total_strain_at_quadrature_points()
#        return_map_smooth_cone_or_apex()
#        assemble_residual_and_consistent_tangent()
#        solve_for_displacement_correction()
#        update_displacement()
#        exit_if_converged()
#      commit_plastic_strain()
#      calculate_footing_reaction()
#      write_pressure_settlement_result()
#  ------------------------------------------------------------
# ==============================================================================


ELASTIC = np.int8(0)
SMOOTH = np.int8(1)
APEX = np.int8(2)

MANDEL_IDENTITY = np.array([1.0, 1.0, 1.0, 0.0])
VOLUMETRIC_PROJECTOR = np.outer(MANDEL_IDENTITY, MANDEL_IDENTITY)
DEVIATORIC_PROJECTOR = np.eye(4) - VOLUMETRIC_PROJECTOR / 3.0


@dataclass(frozen=True)
class Material:
    """Small-strain elastic/perfectly-plastic material parameters."""

    young: float = 1.0e7
    poisson: float = 0.48
    cohesion: float = 450.0
    friction_angle_deg: float = 20.0

    @property
    def shear(self) -> float:
        return self.young / (2.0 * (1.0 + self.poisson))

    @property
    def bulk(self) -> float:
        return self.young / (3.0 * (1.0 - 2.0 * self.poisson))

    @property
    def eta(self) -> float:
        tangent = np.tan(np.deg2rad(self.friction_angle_deg))
        return 3.0 * tangent / np.sqrt(9.0 + 12.0 * tangent**2)

    @property
    def yield_intercept(self) -> float:
        tangent = np.tan(np.deg2rad(self.friction_angle_deg))
        return 3.0 * self.cohesion / np.sqrt(9.0 + 12.0 * tangent**2)

    @property
    def elastic_tangent(self) -> np.ndarray:
        return (
            2.0 * self.shear * DEVIATORIC_PROJECTOR
            + self.bulk * VOLUMETRIC_PROJECTOR
        )


@dataclass
class ConstitutiveResult:
    stress: np.ndarray
    tangent: np.ndarray
    plastic_strain: np.ndarray
    plastic_multiplier: np.ndarray
    branch: np.ndarray


def drucker_prager_yield(stress: np.ndarray, material: Material) -> np.ndarray:
    """Evaluate ``sqrt(J2) + eta*p - c`` for Mandel stress vectors."""

    stress = np.asarray(stress, dtype=float)
    deviatoric_stress = np.einsum(
        "ij,...j->...i", DEVIATORIC_PROJECTOR, stress
    )
    sqrt_j2 = np.linalg.norm(deviatoric_stress, axis=-1) / np.sqrt(2.0)
    mean_stress = np.einsum("i,...i->...", MANDEL_IDENTITY, stress) / 3.0
    return sqrt_j2 + material.eta * mean_stress - material.yield_intercept


def integrate_drucker_prager(
    total_strain: np.ndarray,
    plastic_strain_old: np.ndarray,
    material: Material,
    yield_tolerance: float = 1.0e-12,
) -> ConstitutiveResult:
    """Backward-Euler return mapping at one or more integration points.

    Parameters are arrays whose final axis contains the four plane-strain Mandel
    components.  The returned tangent has shape ``(..., 4, 4)``.  Branch values
    are 0 (elastic), 1 (smooth cone), and 2 (apex).
    """

    strain = np.asarray(total_strain, dtype=float)
    plastic_old = np.asarray(plastic_strain_old, dtype=float)
    if strain.shape != plastic_old.shape or strain.shape[-1] != 4:
        raise ValueError(
            "total_strain and plastic_strain_old must have identical (..., 4) shapes"
        )

    original_shape = strain.shape[:-1]
    strain_2d = strain.reshape((-1, 4))
    plastic_old_2d = plastic_old.reshape((-1, 4))

    shear = material.shear
    bulk = material.bulk
    eta = material.eta
    intercept = material.yield_intercept
    elastic_tangent = material.elastic_tangent

    trial_elastic_strain = strain_2d - plastic_old_2d
    trial_deviatoric_strain = trial_elastic_strain @ DEVIATORIC_PROJECTOR
    deviatoric_norm = np.sqrt(
        np.maximum(
            0.0,
            np.einsum(
                "ni,ni->n", trial_elastic_strain, trial_deviatoric_strain
            ),
        )
    )
    rho_trial = 2.0 * shear * deviatoric_norm
    pressure_trial = bulk * (trial_elastic_strain @ MANDEL_IDENTITY)
    stress_trial = trial_elastic_strain @ elastic_tangent

    apex_denominator = bulk * eta**2
    smooth_denominator = shear + apex_denominator
    criterion_1 = rho_trial / np.sqrt(2.0) + eta * pressure_trial - intercept
    criterion_2 = (
        eta * pressure_trial
        - apex_denominator * rho_trial / (shear * np.sqrt(2.0))
        - intercept
    )

    plastic = criterion_1 > yield_tolerance
    smooth = plastic & (criterion_2 <= 0.0)
    apex = plastic & (criterion_2 > 0.0)

    stress = stress_trial.copy()
    tangent = np.broadcast_to(elastic_tangent, (strain_2d.shape[0], 4, 4)).copy()
    plastic_strain = plastic_old_2d.copy()
    plastic_multiplier = np.zeros(strain_2d.shape[0])
    branch = np.full(strain_2d.shape[0], ELASTIC, dtype=np.int8)

    if np.any(smooth):
        indices = np.flatnonzero(smooth)
        multiplier = criterion_1[indices] / smooth_denominator
        normal = (
            trial_deviatoric_strain[indices]
            / deviatoric_norm[indices, np.newaxis]
        )
        correction_direction = (
            np.sqrt(2.0) * shear * normal
            + bulk * eta * MANDEL_IDENTITY[np.newaxis, :]
        )

        stress[indices] -= multiplier[:, np.newaxis] * correction_direction
        plastic_strain[indices] += multiplier[:, np.newaxis] * (
            normal / np.sqrt(2.0)
            + eta * MANDEL_IDENTITY[np.newaxis, :] / 3.0
        )

        curvature_factor = (
            2.0
            * np.sqrt(2.0)
            * shear**2
            * multiplier
            / rho_trial[indices]
        )
        tangent[indices] -= curvature_factor[:, np.newaxis, np.newaxis] * (
            DEVIATORIC_PROJECTOR[np.newaxis, :, :]
            - np.einsum("ni,nj->nij", normal, normal)
        )
        tangent[indices] -= (
            np.einsum("ni,nj->nij", correction_direction, correction_direction)
            / smooth_denominator
        )

        plastic_multiplier[indices] = multiplier
        branch[indices] = SMOOTH

    if np.any(apex):
        if eta <= 0.0:
            raise RuntimeError("An apex return requires a positive friction parameter")
        indices = np.flatnonzero(apex)
        multiplier = (eta * pressure_trial[indices] - intercept) / apex_denominator
        stress[indices] = (
            intercept / eta * MANDEL_IDENTITY[np.newaxis, :]
        )
        tangent[indices] = 0.0
        plastic_strain[indices] = strain_2d[indices] - (
            intercept
            / (3.0 * bulk * eta)
            * MANDEL_IDENTITY[np.newaxis, :]
        )
        plastic_multiplier[indices] = multiplier
        branch[indices] = APEX

    return ConstitutiveResult(
        stress=stress.reshape(original_shape + (4,)),
        tangent=tangent.reshape(original_shape + (4, 4)),
        plastic_strain=plastic_strain.reshape(original_shape + (4,)),
        plastic_multiplier=plastic_multiplier.reshape(original_shape),
        branch=branch.reshape(original_shape),
    )


def run_constitutive_self_tests(verbose: bool = True) -> None:
    """Check elasticity, both return branches, consistency, and the tangent."""

    material = Material()
    zero = np.zeros(4)

    elastic_strain = np.array([-1.0e-6, -1.0e-6, 0.0, 0.0])
    elastic = integrate_drucker_prager(elastic_strain, zero, material)
    np.testing.assert_equal(elastic.branch, ELASTIC)
    np.testing.assert_allclose(
        elastic.stress, material.elastic_tangent @ elastic_strain, rtol=1e-13
    )

    smooth_strain = np.array([0.0, 0.0, 0.0, 5.0e-4])
    smooth = integrate_drucker_prager(smooth_strain, zero, material)
    np.testing.assert_equal(smooth.branch, SMOOTH)
    np.testing.assert_allclose(
        drucker_prager_yield(smooth.stress, material), 0.0, atol=1e-10
    )
    np.testing.assert_allclose(
        smooth.stress,
        material.elastic_tangent @ (smooth_strain - smooth.plastic_strain),
        rtol=1e-12,
        atol=1e-9,
    )

    apex_strain = np.array([5.0e-4, 5.0e-4, 0.0, 0.0])
    apex = integrate_drucker_prager(apex_strain, zero, material)
    np.testing.assert_equal(apex.branch, APEX)
    np.testing.assert_allclose(
        drucker_prager_yield(apex.stress, material), 0.0, atol=1e-10
    )
    np.testing.assert_allclose(
        DEVIATORIC_PROJECTOR @ apex.stress, 0.0, atol=1e-10
    )
    np.testing.assert_allclose(apex.tangent, 0.0, atol=0.0)

    direction = np.array([0.17, -0.31, 0.09, 0.41])
    step = 1.0e-8
    stress_plus = integrate_drucker_prager(
        smooth_strain + step * direction, zero, material
    ).stress
    stress_minus = integrate_drucker_prager(
        smooth_strain - step * direction, zero, material
    ).stress
    finite_difference = (stress_plus - stress_minus) / (2.0 * step)
    tangent_action = smooth.tangent @ direction
    relative_error = np.linalg.norm(
        finite_difference - tangent_action
    ) / np.linalg.norm(tangent_action)
    if relative_error > 2.0e-7:
        raise AssertionError(
            f"Consistent-tangent finite-difference error is {relative_error:.3e}"
        )

    if verbose:
        print("Constitutive checks passed")
        print(f"  eta               = {material.eta:.12g}")
        print(f"  yield intercept   = {material.yield_intercept:.12g}")
        smooth_yield = drucker_prager_yield(smooth.stress, material)
        print(f"  smooth yield      = {smooth_yield:.3e}")
        print(f"  tangent rel. error= {relative_error:.3e}")


def run_strip_footing(
    material: Material,
    cells_per_side: int = 20,
    cell_type: str = "triangle",
    mesh_file: Path | None = None,
    element_family: str = "auto",
    degree: int = 2,
    quadrature_degree: int | None = None,
    increments: int = 12,
    maximum_settlement: float = 0.03,
    maximum_newton_iterations: int = 30,
    relative_tolerance: float = 1.0e-8,
    paraview_output: Path | None = None,
) -> list[dict[str, float]]:
    """Solve the displacement-controlled plane-strain strip-footing problem."""

    import basix
    import basix.ufl
    import ufl
    from dolfinx import fem, io, mesh
    import dolfinx.fem.petsc
    from mpi4py import MPI
    from petsc4py import PETSc

    comm = MPI.COMM_WORLD
    length = 10.0
    footing_width = 1.0
    if mesh_file is None:
        dolfinx_cell_type = {
            "triangle": mesh.CellType.triangle,
            "quadrilateral": mesh.CellType.quadrilateral,
        }[cell_type]
        domain = mesh.create_rectangle(
            comm,
            [np.array([0.0, 0.0]), np.array([length, length])],
            [cells_per_side, cells_per_side],
            cell_type=dolfinx_cell_type,
        )
    else:
        domain = io.gmsh.read_from_msh(mesh_file, comm, rank=0, gdim=2).mesh

    is_triangle = domain.topology.cell_type == mesh.CellType.triangle
    if element_family == "auto":
        element_family = "Lagrange" if is_triangle else "S"
    if is_triangle and element_family == "S":
        raise ValueError("Serendipity elements are not defined on triangles")
    if quadrature_degree is None:
        # Seven-point triangle and 3x3 quadrilateral rules from the published
        # P2 and Q2 implementations, respectively.
        quadrature_degree = 5 if is_triangle else 4
    vector_space = fem.functionspace(domain, (element_family, degree, (2,)))

    def bottom(x: np.ndarray) -> np.ndarray:
        return np.isclose(x[1], 0.0)

    def left(x: np.ndarray) -> np.ndarray:
        return np.isclose(x[0], 0.0)

    def right(x: np.ndarray) -> np.ndarray:
        return np.isclose(x[0], length)

    def footing(x: np.ndarray) -> np.ndarray:
        return np.isclose(x[1], length) & (x[0] <= footing_width + 1.0e-12)

    scalar_x_space, _ = vector_space.sub(0).collapse()
    scalar_y_space, _ = vector_space.sub(1).collapse()
    facet_dimension = domain.topology.dim - 1
    bottom_facets = mesh.locate_entities_boundary(domain, facet_dimension, bottom)
    left_facets = mesh.locate_entities_boundary(domain, facet_dimension, left)
    right_facets = mesh.locate_entities_boundary(domain, facet_dimension, right)
    footing_facets = mesh.locate_entities_boundary(domain, facet_dimension, footing)
    bottom_y_dofs = fem.locate_dofs_topological(
        (vector_space.sub(1), scalar_y_space), facet_dimension, bottom_facets
    )
    left_x_dofs = fem.locate_dofs_topological(
        (vector_space.sub(0), scalar_x_space), facet_dimension, left_facets
    )
    right_x_dofs = fem.locate_dofs_topological(
        (vector_space.sub(0), scalar_x_space), facet_dimension, right_facets
    )
    footing_y_dofs = fem.locate_dofs_topological(
        (vector_space.sub(1), scalar_y_space), facet_dimension, footing_facets
    )

    zero_x = fem.Function(scalar_x_space)
    zero_y = fem.Function(scalar_y_space)
    prescribed_vertical_displacement = fem.Function(scalar_y_space)
    boundary_conditions = [
        fem.dirichletbc(zero_y, bottom_y_dofs, vector_space.sub(1)),
        fem.dirichletbc(zero_x, left_x_dofs, vector_space.sub(0)),
        fem.dirichletbc(zero_x, right_x_dofs, vector_space.sub(0)),
        fem.dirichletbc(
            prescribed_vertical_displacement,
            footing_y_dofs,
            vector_space.sub(1),
        ),
    ]

    quadrature_element_vector = basix.ufl.quadrature_element(
        domain.basix_cell(),
        value_shape=(4,),
        scheme="default",
        degree=quadrature_degree,
    )
    quadrature_element_tensor = basix.ufl.quadrature_element(
        domain.basix_cell(),
        value_shape=(4, 4),
        scheme="default",
        degree=quadrature_degree,
    )
    quadrature_element_scalar = basix.ufl.quadrature_element(
        domain.basix_cell(),
        value_shape=(),
        scheme="default",
        degree=quadrature_degree,
    )
    state_space = fem.functionspace(domain, quadrature_element_vector)
    tangent_space = fem.functionspace(domain, quadrature_element_tensor)
    scalar_state_space = fem.functionspace(domain, quadrature_element_scalar)

    stress = fem.Function(state_space, name="Mandel_stress")
    plastic_strain_old = fem.Function(state_space, name="Plastic_strain_old")
    plastic_strain_candidate = fem.Function(state_space, name="Plastic_strain")
    consistent_tangent = fem.Function(tangent_space, name="Consistent_tangent")
    plastic_multiplier = fem.Function(scalar_state_space, name="Plastic_multiplier")
    return_branch = fem.Function(scalar_state_space, name="Return_branch")

    displacement = fem.Function(vector_space, name="Displacement")
    correction = fem.Function(vector_space, name="Newton_correction")
    trial = ufl.TrialFunction(vector_space)
    test = ufl.TestFunction(vector_space)

    def strain_mandel(vector: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        strain = ufl.sym(ufl.grad(vector))
        return ufl.as_vector(
            [strain[0, 0], strain[1, 1], 0.0, np.sqrt(2.0) * strain[0, 1]]
        )

    metadata = {
        "quadrature_degree": quadrature_degree,
        "quadrature_scheme": "default",
    }
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)
    residual = ufl.dot(strain_mandel(test), stress) * dx
    jacobian = ufl.dot(
        strain_mandel(test), ufl.dot(consistent_tangent, strain_mandel(trial))
    ) * dx
    residual_form = fem.form(residual)
    negative_residual_form = fem.form(-residual)
    jacobian_form = fem.form(jacobian)

    basix_cell_type = getattr(basix.CellType, domain.topology.cell_type.name)
    quadrature_points, _ = basix.make_quadrature(
        basix_cell_type, quadrature_degree
    )
    cell_map = domain.topology.index_map(domain.topology.dim)
    cells = np.arange(cell_map.size_local + cell_map.num_ghosts, dtype=np.int32)
    strain_expression = fem.Expression(strain_mandel(displacement), quadrature_points)

    def update_constitutive_state() -> ConstitutiveResult:
        strain_values = strain_expression.eval(domain, cells).reshape((-1, 4))
        old_values = plastic_strain_old.x.array.reshape((-1, 4))
        if strain_values.shape != old_values.shape:
            raise RuntimeError(
                "Quadrature expression and state-function layouts do not agree: "
                f"{strain_values.shape} != {old_values.shape}"
            )
        result = integrate_drucker_prager(strain_values, old_values, material)
        stress.x.array[:] = result.stress.reshape(-1)
        plastic_strain_candidate.x.array[:] = result.plastic_strain.reshape(-1)
        consistent_tangent.x.array[:] = result.tangent.reshape(-1)
        plastic_multiplier.x.array[:] = result.plastic_multiplier.reshape(-1)
        return_branch.x.array[:] = result.branch.reshape(-1)
        return result

    solver = PETSc.KSP().create(comm)
    solver.setType("preonly")
    solver.getPC().setType("lu")

    paraview_writer = None
    if paraview_output is not None:
        if comm.rank == 0:
            paraview_output.parent.mkdir(parents=True, exist_ok=True)
        comm.barrier()
        # VTX/BP preserves the quadratic displacement field and is read directly
        # by ParaView's ADIOS2 reader.
        paraview_writer = io.VTXWriter(comm, paraview_output, displacement)
        paraview_writer.write(0.0)

    results: list[dict[str, float]] = []
    settlements = np.linspace(0.0, maximum_settlement, increments + 1)[1:]
    for step_number, settlement in enumerate(settlements, start=1):
        prescribed_vertical_displacement.x.array[:] = PETSc.ScalarType(-settlement)
        prescribed_vertical_displacement.x.scatter_forward()
        initial_residual_norm = None
        converged = False

        for iteration in range(maximum_newton_iterations + 1):
            constitutive_result = update_constitutive_state()

            rhs = dolfinx.fem.petsc.assemble_vector(negative_residual_form)
            dolfinx.fem.petsc.apply_lifting(
                rhs,
                [jacobian_form],
                bcs=[boundary_conditions],
                x0=[displacement.x.petsc_vec],
                alpha=1.0,
            )
            rhs.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES,
                mode=PETSc.ScatterMode.REVERSE,
            )
            dolfinx.fem.petsc.set_bc(
                rhs,
                boundary_conditions,
                x0=displacement.x.petsc_vec,
                alpha=1.0,
            )
            residual_norm = rhs.norm()
            if initial_residual_norm is None:
                initial_residual_norm = max(residual_norm, 1.0)

            if iteration > 0 and residual_norm <= max(
                1.0e-10, relative_tolerance * initial_residual_norm
            ):
                converged = True
                rhs.destroy()
                break
            if iteration == maximum_newton_iterations:
                rhs.destroy()
                break

            matrix = dolfinx.fem.petsc.assemble_matrix(
                jacobian_form, bcs=boundary_conditions
            )
            matrix.assemble()
            correction.x.array[:] = 0.0
            solver.setOperators(matrix)
            solver.solve(rhs, correction.x.petsc_vec)
            correction.x.scatter_forward()
            displacement.x.petsc_vec.axpy(1.0, correction.x.petsc_vec)
            displacement.x.scatter_forward()
            matrix.destroy()
            rhs.destroy()

        if not converged:
            raise RuntimeError(
                f"Newton solver failed at settlement {settlement:.6g}; "
                f"last residual {residual_norm:.6e}"
            )

        # Commit history only after the global equilibrium iteration converges.
        plastic_strain_old.x.array[:] = plastic_strain_candidate.x.array

        internal_force = dolfinx.fem.petsc.assemble_vector(residual_form)
        internal_force.ghostUpdate(
            addv=PETSc.InsertMode.ADD_VALUES,
            mode=PETSc.ScatterMode.REVERSE,
        )
        owned_size = (
            vector_space.dofmap.index_map.size_local
            * vector_space.dofmap.index_map_bs
        )
        # For a (parent subspace, collapsed space) pair, DOLFINx returns the
        # corresponding parent and collapsed dof arrays as a two-item list.
        footing_parent_dofs = np.asarray(footing_y_dofs[0])
        owned_footing_dofs = footing_parent_dofs[footing_parent_dofs < owned_size]
        local_vertical_reaction = np.sum(internal_force.array[owned_footing_dofs])
        vertical_reaction = comm.allreduce(local_vertical_reaction, op=MPI.SUM)
        normalized_pressure = -vertical_reaction / (
            footing_width * material.cohesion
        )
        owned_quadrature_points = cell_map.size_local * len(quadrature_points)
        owned_branches = constitutive_result.branch[:owned_quadrature_points]
        local_smooth = np.count_nonzero(owned_branches == SMOOTH)
        local_apex = np.count_nonzero(owned_branches == APEX)
        smooth_count = comm.allreduce(local_smooth, op=MPI.SUM)
        apex_count = comm.allreduce(local_apex, op=MPI.SUM)

        row = {
            "step": float(step_number),
            "settlement": float(settlement),
            "normalized_pressure": float(normalized_pressure),
            "newton_iterations": float(iteration),
            "relative_residual": float(residual_norm / initial_residual_norm),
            "smooth_points": float(smooth_count),
            "apex_points": float(apex_count),
        }
        results.append(row)
        if comm.rank == 0:
            print(
                f"step={step_number:3d} settlement={settlement:.6e} "
                f"pressure/c0={normalized_pressure:.8e} Newton={iteration:2d} "
                f"relres={row['relative_residual']:.3e} "
                f"smooth={smooth_count} apex={apex_count}"
            )
        if paraview_writer is not None:
            paraview_writer.write(float(settlement))
        internal_force.destroy()

    if paraview_writer is not None:
        paraview_writer.close()
    solver.destroy()
    return results


def write_results(path: Path, results: list[dict[str, float]]) -> None:
    if not results:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=list(results[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-test-only",
        action="store_true",
        help="run material-point checks without constructing a DOLFINx problem",
    )
    parser.add_argument(
        "--cells",
        type=int,
        default=20,
        help="cells per side (use a multiple of 10 to align the footing edge)",
    )
    parser.add_argument(
        "--cell-type",
        choices=("triangle", "quadrilateral"),
        default="triangle",
        help="cell type used when --mesh is not supplied",
    )
    parser.add_argument(
        "--mesh",
        type=Path,
        help="optional Gmsh mesh; overrides --cells and --cell-type",
    )
    parser.add_argument("--degree", type=int, default=2, help="displacement degree")
    parser.add_argument(
        "--element-family",
        choices=("auto", "S", "Lagrange"),
        default="auto",
        help="auto selects P2 on triangles and eight-node Q2 on quadrilaterals",
    )
    parser.add_argument(
        "--quadrature-degree",
        type=int,
        help="default: degree 5 on triangles and degree 4 on quadrilaterals",
    )
    parser.add_argument("--increments", type=int, default=12)
    parser.add_argument("--maximum-settlement", type=float, default=0.03)
    parser.add_argument("--maximum-newton-iterations", type=int, default=30)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--paraview-output",
        type=Path,
        help="write the quadratic displacement time series as a ParaView VTX/BP file",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    try:
        from mpi4py import MPI

        is_root = MPI.COMM_WORLD.rank == 0
    except ImportError:
        is_root = True
    run_constitutive_self_tests(verbose=is_root)
    if arguments.self_test_only:
        return

    results = run_strip_footing(
        Material(),
        cells_per_side=arguments.cells,
        cell_type=arguments.cell_type,
        mesh_file=arguments.mesh,
        element_family=arguments.element_family,
        degree=arguments.degree,
        quadrature_degree=arguments.quadrature_degree,
        increments=arguments.increments,
        maximum_settlement=arguments.maximum_settlement,
        maximum_newton_iterations=arguments.maximum_newton_iterations,
        relative_tolerance=arguments.relative_tolerance,
        paraview_output=arguments.paraview_output,
    )
    if arguments.output is not None:
        if is_root:
            write_results(arguments.output, results)
            print(f"Wrote {arguments.output}")


if __name__ == "__main__":
    main()
