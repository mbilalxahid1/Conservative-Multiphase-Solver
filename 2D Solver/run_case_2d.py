# Main 2D driver. No 1D run option is present in this file.

import glob
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from case_2D_ShockBubbleCollapse import get_config
from models.mesh_2d import build_mesh_2d
from models.state_vars import allocate_state
from models.eos_nasg import primitive_to_conservative_2d
from models.geometry_2d import get_initial_geometry, initial_state_at_point_2d
from models.time_step_cfl import compute_dt_cfl
from models.rhs_euler_2d import rhs_euler_2d
from models.level_set_solution import (
    initialize_phi_and_material_id_2d,
    level_set_solution,
)
from models.bi_material_cell_updates import update_bi_material_cells_2d
from models.io_outputs import write_snapshot
from models.diagnostic_outputs_2d import DiagnosticOutput2D
from models.parallel_numba import configure_numba, numba_enabled_for_size


def _snapshot_extras(output_cfg, material_id, phi):
    extras = {}
    fields = output_cfg.get("fields", {})
    if fields.get("material_id", False):
        extras["material_id"] = material_id
    if fields.get("levelset0", False):
        extras["levelset0"] = phi
    return extras


def _validate_solver_options(cfg):
    mesh_type = str(cfg["mesh"].get("type", "TwoDimensional")).strip().lower()
    if mesh_type not in ("twodimensional", "2d", "two_dimensional"):
        raise ValueError("run_case_2d supports 2D only.")

    method = str(
        cfg.get("multiphase", {}).get("phase_change", "Extrapolation")
    ).strip().lower()
    if method not in ("volume", "extrapolation"):
        raise ValueError(
            "run_case_2d: phase_change must be Extrapolation or Volume.")

    flux = str(
        cfg.get("space", {}).get("navier_stokes", {}).get("flux", "Rusanov")
    ).strip().lower()
    if flux not in ("rusanov", "roe", "llf"):
        raise ValueError("run_case_2d: flux must be Rusanov or Roe.")
    from models.rhs_euler_2d import _reconstruction_type
    reconstruction = _reconstruction_type(
        cfg.get("space", {}).get("navier_stokes", {}).get("reconstruction"))
    if reconstruction not in ("constant", "linear"):
        raise ValueError("run_case_2d: reconstruction must be Constant or Linear.")

    if method == "volume":
        from models.level_set_solution import _read_level_set_settings
        background, interface, _ = _read_level_set_settings(cfg)
        if (background, interface) != (0, 1) or sorted(cfg["equations"]["materials"]) != [0, 1]:
            raise NotImplementedError("2D Volume currently requires exactly material IDs 0 and 1.")

    time_type = str(cfg.get("time", {}).get("type", "Explicit")).strip().lower()
    rk_type = str(
        cfg.get("time", {}).get("explicit", {}).get("type", "RungeKutta2")
    ).strip().lower()
    if time_type != "explicit" or rk_type not in ("rungekutta2", "rk2"):
        raise NotImplementedError(
            "run_case_2d: this 2D stage uses explicit SSP-RK2 only.")


def _run_volume(cfg, U_initial, phi, material_id, x, x_faces, dx,
                y, y_faces, dy, parallel_info, max_steps, quiet):
    """Regular simulation time/output loop using authoritative sub-volume Q."""
    from models.bi_material_state_calculator.volume import (
        VolumeStateError, initialize_volume_state_2d,
        reconstruct_volume_geometry_2d, recover_volume_states_2d,
        representative_volume_state_2d, compute_volume_dt_cfl_2d,
        volume_step_2d,
    )
    materials = cfg["equations"]["materials"]
    geometry = reconstruct_volume_geometry_2d(phi, x_faces, y_faces, cfg)
    Q0, Q1 = initialize_volume_state_2d(
        U_initial, material_id, geometry, materials, x, y)

    output_cfg = dict(cfg["output"])
    folder = str(output_cfg.get("folder_name", output_cfg.get("prefix", "results_2d/")))
    output_cfg["folder_name"] = folder
    output_cfg["prefix"] = folder
    os.makedirs(folder, exist_ok=True)
    solution_name = str(output_cfg.get("solution_name", "solution"))
    for pattern in (f"{solution_name}_*.vtr", f"{solution_name}.pvd"):
        for file_name in glob.glob(os.path.join(folder, pattern)):
            os.remove(file_name)

    diagnostics = DiagnosticOutput2D(output_cfg, x, y, dx, dy, materials)
    t = 0.0
    step = 0
    final_time = float(cfg["time"]["max_time"])
    CFL = float(cfg["time"]["cfl"])
    output_interval = float(output_cfg.get("time_interval", final_time))
    next_output = output_interval if output_interval > 0 else np.inf
    options = cfg.get("multiphase", {}).get("volume_solver", {})
    max_retries = int(options.get("max_retries", 8))
    if max_retries < 0:
        raise ValueError("Volume max_retries must be nonnegative.")
    adaptive_enabled = bool(options.get("adaptive_merge_enabled", False))
    adaptive_after = int(options.get("adaptive_merge_after_retries", 3))
    adaptive_max_cells = int(options.get("adaptive_merge_max_cells", 5))
    if adaptive_enabled and (adaptive_after < 0 or not 2 <= adaptive_max_cells <= 5):
        raise ValueError("Adaptive merge requires after_retries>=0 and 2<=max_cells<=5.")
    pre_riemann_enabled = bool(options.get("pre_riemann_merge_enabled", False))
    pre_riemann_after = int(options.get("pre_riemann_merge_after_retries", 3))
    pre_riemann_max_cells = int(options.get("pre_riemann_merge_max_cells", 5))
    if pre_riemann_enabled and (pre_riemann_after < 0 or
                                not 2 <= pre_riemann_max_cells <= 5):
        raise ValueError("Pre-Riemann merge requires after_retries>=0 and 2<=max_cells<=5.")
    extrapolation_enabled = bool(options.get("conservative_extrapolation_enabled", False))
    extrapolation_after = int(options.get("conservative_extrapolation_after_retries", 3))
    if extrapolation_enabled and extrapolation_after < 0:
        raise ValueError("Conservative extrapolation requires after_retries>=0.")

    def display_fields(q0, q1, geom, phi_local):
        u0, u1 = recover_volume_states_2d(
            q0, q1, geom, materials, parallel_cfg=cfg)
        ids = (phi_local < 0).astype(int)
        # A cell-center ID is only meaningful when that material has area.
        # Otherwise use the reconstructed majority to avoid invalid U=0.
        missing = ((ids == 0) & (geom["area_0"] <= 0)
                   | (ids == 1) & (geom["area_1"] <= 0))
        ids[missing] = (geom["area_1"][missing] > geom["area_0"][missing])
        u = representative_volume_state_2d(u0, u1, geom, material_id=ids)
        return u, ids

    U, display_ids = display_fields(Q0, Q1, geometry, phi)
    volume_state = (Q0, Q1, geometry)
    write_snapshot(output_cfg, 0, t, x_faces, y_faces, U,
                   _snapshot_extras(output_cfg, display_ids, phi),
                   materials, display_ids, volume_state=volume_state)
    diagnostics.initialize(0, t, U, display_ids, phi, volume_state=volume_state)
    history = [diagnostics.record_volume_budget(0, t, 0.0, Q0, Q1, geometry)]

    while t < final_time and step < max_steps:
        dt = compute_volume_dt_cfl_2d(
            Q0, Q1, geometry, dx, dy, materials, CFL, parallel_cfg=cfg)
        dt = min(dt, final_time - t)
        retries = 0
        while True:
            # Three attempts use the original two-cell merge: dt0, dt0/2,
            # and dt0/4. After three failed attempts, the third halving
            # starts adaptive merging at dt0/8. Later failures still halve dt.
            adaptive_now = adaptive_enabled and retries >= adaptive_after
            extrapolation_now = extrapolation_enabled and retries >= extrapolation_after
            pre_riemann_now = pre_riemann_enabled and retries >= pre_riemann_after
            try:
                step_args = {}
                if adaptive_now:
                    step_args["adaptive_merge_max_cells"] = adaptive_max_cells
                if extrapolation_now:
                    step_args["conservative_extrapolation_enabled"] = True
                if pre_riemann_now:
                    step_args["pre_riemann_merge_enabled"] = True
                    step_args["pre_riemann_merge_max_cells"] = pre_riemann_max_cells
                q0_new, q1_new, phi_new, ids_new, info = volume_step_2d(
                    Q0, Q1, phi, x_faces, y_faces, cfg, dt,
                    step=step + 1, **step_args)
                break
            except VolumeStateError as exc:
                if retries >= max_retries or t + 0.5*dt <= t:
                    raise RuntimeError(
                        f"2D Volume step {step+1} failed after {retries} retries "
                        f"at time {t:.8e} (dt={dt:.4e}): {exc}") from exc
                retries += 1
                dt *= 0.5
                if not quiet:
                    print(f"Volume retry: step={step+1}, dt={dt:.4e}, reason={exc}")
                    if adaptive_enabled and retries == adaptive_after:
                        print(f"Volume adaptive merge activated: step={step+1}, "
                              f"max_cells={adaptive_max_cells}, retry={retries}")
                    if extrapolation_enabled and retries == extrapolation_after:
                        print(f"Volume conservative extrapolation activated: "
                              f"step={step+1}, retry={retries}")
                    if pre_riemann_enabled and retries == pre_riemann_after:
                        print(f"Volume pre-Riemann conservative merge activated: "
                              f"step={step+1}, max_cells={pre_riemann_max_cells}, "
                              f"retry={retries}")

        U_old = U
        display_ids_old = display_ids
        Q0, Q1, phi, material_id = q0_new, q1_new, phi_new, ids_new
        geometry = reconstruct_volume_geometry_2d(phi, x_faces, y_faces, cfg)
        U, display_ids = display_fields(Q0, Q1, geometry, phi)
        t += dt
        step += 1
        crossed = display_ids != display_ids_old
        volume_state = (Q0, Q1, geometry)
        # Add this step's boundary flux before writing the integrations.
        # record_volume_budget also maintains volume_conservation.txt.
        history.append(diagnostics.record_volume_budget(
            step, t, dt, Q0, Q1, geometry, info))
        diagnostics.record_step(
            step, t, dt, U_old, U, U,
            display_ids_old, display_ids, crossed,
            phi=phi, reinitialized=info["reinitialized"],
            volume_state=volume_state)
        if not quiet:
            print(f"Volume step {step}: t={t:.8e}, dt={dt:.3e}, "
                  f"retries={retries}, mass_error="
                  f"{max(abs(history[-1]['mass_error'])):.3e}")
            pre_events = (info.get("pre_riemann_merge_events_rk1", []) +
                          info.get("pre_riemann_merge_events_rk2", []))
            if pre_events:
                print(f"Volume pre-Riemann merge used: step={step}, "
                      f"groups={len(pre_events)}, "
                      f"largest_group={max(event['group_cells'] for event in pre_events)}, "
                      f"segments={[event['segment'] for event in pre_events]}, "
                      f"first_cells_(material,i,j)="
                      f"{[(event['material'], *event['cell']) for event in pre_events[:8]]}")
            event_reports = (info["rk1_events"], info["rk2_events"],
                             info["reinit_events"])
            adaptive_groups = sum(int(item.get("adaptive_merge_groups", 0))
                                  for item in event_reports)
            if adaptive_groups:
                max_group = max(int(item.get("adaptive_merge_max_group_cells", 0))
                                for item in event_reports)
                print(f"Volume adaptive merge used: step={step}, "
                      f"groups={adaptive_groups}, largest_group={max_group} cells")
            extrap_repaired = sum(int(item.get("extrapolation_repaired", 0))
                                  for item in event_reports)
            if extrap_repaired:
                donors = max(int(item.get("extrapolation_max_donors_used", 0))
                             for item in event_reports)
                cells = [pos for item in event_reports
                         for pos in item.get("extrapolation_cells", [])]
                min_blend = min(float(item.get("extrapolation_min_blend", 1.0))
                                for item in event_reports)
                print(f"Volume conservative extrapolation used: step={step}, "
                      f"repaired={extrap_repaired}, largest_donor_group={donors} donors, "
                      f"min_blend={min_blend:.6f}, "
                      f"first_cells_(material,i,j,donors)={cells[:8]}")
            if info.get("orphan_discarded_groups", 0):
                loss = info["orphan_discarded_Q"]
                absolute = info["orphan_discarded_abs_Q"]
                print(f"Volume orphan cleanup (non-conservative): step={step}, "
                      f"groups={info['orphan_discarded_groups']}, "
                      f"cells={info['orphan_discarded_cells']}, "
                      f"effective_lost_mass_0={loss[0,0]:.16e}, "
                      f"effective_lost_mass_1={loss[1,0]:.16e}, "
                      f"absolute_mass_removed={absolute[:,0].sum():.16e}; "
                      f"details in volume_orphan_discards.txt")
        if t >= next_output - 1.e-14 or t >= final_time - 1.e-14:
            write_snapshot(
                output_cfg, step, t, x_faces, y_faces, U,
                _snapshot_extras(output_cfg, display_ids, phi),
                materials, display_ids, volume_state=volume_state)
            while next_output <= t + 1.e-14:
                next_output += output_interval

    if step >= max_steps and t < final_time:
        raise RuntimeError("run_case_2d: maximum step count reached before final time.")
    diagnostics.finalize(step, t, U, display_ids, phi, volume_state=volume_state)
    return {
        "U": U, "Q0": Q0, "Q1": Q1, "geometry": geometry,
        "phi": phi, "material_id": material_id,
        "display_material_id": display_ids,
        "x_cells": x, "y_cells": y,
        "x_faces": x_faces, "y_faces": y_faces,
        "dx": dx, "dy": dy, "t": t, "step": step,
        "volume_conservation": history,
        "rk_mismatch_x": np.zeros((len(y), len(x) + 1, 4)),
        "rk_mismatch_y": np.zeros((len(y) + 1, len(x), 4)),
        "parallel": parallel_info,
    }


def run(cfg=None, max_steps=1000000, quiet=False):
    cfg = get_config() if cfg is None else cfg
    _validate_solver_options(cfg)

    parallel_info = configure_numba(cfg)
    parallel_cfg = cfg.get("parallel", {})
    if not quiet:
        if parallel_info["active"]:
            print(
                "Numba parallel backend: enabled, "
                f"threads={parallel_info['num_threads']}, "
                f"min_cells={parallel_info['min_cells']}")
        else:
            print("Numba parallel backend: disabled")

    x, x_faces, dx, y, y_faces, dy = build_mesh_2d(cfg["mesh"])
    nx = len(x)
    ny = len(y)
    parallel_info["used_for_mesh"] = numba_enabled_for_size(
        parallel_cfg, nx * ny)
    volume_mode = str(cfg.get("multiphase", {}).get(
        "phase_change", "Extrapolation")).strip().lower() == "volume"
    if volume_mode:
        # Volume has dedicated parallel face, Riemann, state and level-set
        # kernels, using the same activation threshold as Extrapolation.
        parallel_info["volume_numba"] = parallel_info["used_for_mesh"]
    if not quiet:
        print(f"2D mesh: nx={nx}, ny={ny}, cells={nx * ny}")
        if parallel_info["active"]:
            print(
                "Numba hot kernels active for this mesh: "
                f"{parallel_info['used_for_mesh']}")

    materials = cfg["equations"]["materials"]
    phi, material_id = initialize_phi_and_material_id_2d(
        x, y, dx, dy, cfg)

    U = allocate_state(nx, ny)
    default_state = cfg["initial_condition"]["default_state"]
    geometry = get_initial_geometry(cfg)

    # Keep the old rectangular initial_vapor_region path for the Rallu
    # regression case. New Circle/Square cases use geometry_2d.py.
    region = cfg["initial_condition"].get("initial_vapor_region", {})
    region_state = (
        region.get("state", {})
        if isinstance(region, dict) and region.get("enabled", False)
        else {})
    region_material_id = int(region_state.get("material_id", -999999))

    for j in range(ny):
        for i in range(nx):
            if geometry is not None:
                state = initial_state_at_point_2d(x[i], y[j], cfg)
            elif int(material_id[j, i]) == region_material_id:
                state = region_state
            else:
                state = default_state

            mid = int(state["material_id"])
            if mid not in materials:
                raise ValueError(
                    f"run_case_2d: material ID {mid} in the initial condition is undefined.")

            U[j, i] = primitive_to_conservative_2d(
                float(state["rho"]),
                float(state["u"]),
                float(state.get("v", 0.0)),
                float(state["p"]),
                materials[mid])

    if volume_mode:
        return _run_volume(
            cfg, U, phi, material_id, x, x_faces, dx,
            y, y_faces, dy, parallel_info, max_steps, quiet)

    output_cfg = cfg["output"]
    folder = str(output_cfg.get("folder_name", output_cfg.get("prefix", "results_2d/")))
    os.makedirs(folder, exist_ok=True)

    # Clean only the VTK solution series here. DiagnosticOutput2D separately
    # cleans its own text outputs.
    solution_name = str(output_cfg.get("solution_name", "solution"))
    for pattern in (
            f"{solution_name}_*.vtr",
            f"{solution_name}.pvd"):
        for file_name in glob.glob(os.path.join(folder, pattern)):
            try:
                os.remove(file_name)
            except OSError:
                pass

    write_snapshot(
        output_cfg, 0, 0.0, x_faces, y_faces, U,
        _snapshot_extras(output_cfg, material_id, phi),
        materials, material_id)

    diagnostics = DiagnosticOutput2D(
        output_cfg, x, y, dx, dy, materials)
    diagnostics.initialize(0, 0.0, U, material_id, phi)

    t = 0.0
    step = 0
    final_time = float(cfg["time"]["max_time"])
    CFL = float(cfg["time"]["cfl"])
    output_interval = float(output_cfg.get("time_interval", final_time))
    next_output = output_interval if output_interval > 0.0 else np.inf

    reconstruction_cfg = cfg["space"]["navier_stokes"].get(
        "reconstruction", {"type": "Constant"})

    # Per-stage extensive FIVER mismatches. For SSP-RK2, the integrated
    # mismatch over one physical time step is
    #   0.5*dt*loss_RK1 + 0.5*dt*loss_RK2,
    # matching M2C's RK2 bookkeeping.
    mismatch_x_rk1 = np.zeros((ny, nx + 1, 4), dtype=float)
    mismatch_y_rk1 = np.zeros((ny + 1, nx, 4), dtype=float)
    mismatch_x_rk2 = np.zeros_like(mismatch_x_rk1)
    mismatch_y_rk2 = np.zeros_like(mismatch_y_rk1)

    # Reuse the same face-flux buffers for both RK stages. Their boundary
    # values, plus the FIVER mismatches, give the actual interior-side flux.
    face_flux_x = np.empty_like(mismatch_x_rk1)
    face_flux_y = np.empty_like(mismatch_y_rk1)

    last_rk_mismatch_x = np.zeros_like(mismatch_x_rk1)
    last_rk_mismatch_y = np.zeros_like(mismatch_y_rk1)

    while t < final_time and step < max_steps:
        dt = compute_dt_cfl(
            U, dx, dy, materials, material_id, CFL, parallel_cfg)
        dt = min(dt, final_time - t)

        U_old = U.copy()
        material_id_old = material_id.copy()

        callback_rk1 = diagnostics.make_interface_callback(
            step + 1, t, "RK1")
        k1 = rhs_euler_2d(
            U, x_faces, y_faces, dx, dy,
            materials, material_id, reconstruction_cfg, cfg,
            interface_callback=callback_rk1,
            flux_mismatch_x=mismatch_x_rk1,
            flux_mismatch_y=mismatch_y_rk1,
            face_flux_x_output=face_flux_x,
            face_flux_y_output=face_flux_y)
        boundary_rk1 = diagnostics.boundary_flux_from_faces(
            face_flux_x, face_flux_y, mismatch_x_rk1, mismatch_y_rk1,
            material_id_old)
        U_stage = U + dt * k1

        callback_rk2 = diagnostics.make_interface_callback(
            step + 1, t + dt, "RK2")
        k2 = rhs_euler_2d(
            U_stage, x_faces, y_faces, dx, dy,
            materials, material_id, reconstruction_cfg, cfg,
            interface_callback=callback_rk2,
            flux_mismatch_x=mismatch_x_rk2,
            flux_mismatch_y=mismatch_y_rk2,
            face_flux_x_output=face_flux_x,
            face_flux_y_output=face_flux_y)
        boundary_rk2 = diagnostics.boundary_flux_from_faces(
            face_flux_x, face_flux_y, mismatch_x_rk2, mismatch_y_rk2,
            material_id_old)

        U_candidate = 0.5 * U + 0.5 * (U_stage + dt * k2)

        # Extensive state mismatch accumulated by SSP-RK2 during this step.
        last_rk_mismatch_x = 0.5 * dt * (
            mismatch_x_rk1 + mismatch_x_rk2)
        last_rk_mismatch_y = 0.5 * dt * (
            mismatch_y_rk1 + mismatch_y_rk2)

        phi, material_id, reinitialized = level_set_solution(
            U, U_stage, phi, material_id,
            x, y, dx, dy, dt, step + 1, cfg)

        U, crossed = update_bi_material_cells_2d(
            U_candidate, material_id_old, material_id,
            x, y, dx, dy, cfg)

        if not np.all(np.isfinite(U)) or np.any(U[..., 0] <= 0.0):
            raise RuntimeError(
                f"run_case_2d: invalid conservative state after step {step + 1}.")

        t += dt
        step += 1

        # SSP-RK2 time integration of outward mass, momentum and energy flux.
        diagnostics.add_boundary_flux(0.5 * dt * (boundary_rk1 + boundary_rk2))

        diagnostics.record_step(
            step, t, dt,
            U_old, U_candidate, U,
            material_id_old, material_id, crossed,
            phi=phi,
            reinitialized=reinitialized,
            rk_mismatch_x=last_rk_mismatch_x,
            rk_mismatch_y=last_rk_mismatch_y)

        if not quiet and (
                step <= 5 or step % 50 == 0
                or np.any(crossed) or reinitialized):
            print(
                f"Step {step}: t={t:.8e}, dt={dt:.3e}, "
                f"crossed={int(np.count_nonzero(crossed))}, "
                f"reinitialized={reinitialized}")

        if t >= next_output - 1.0e-14 or t >= final_time - 1.0e-14:
            write_snapshot(
                output_cfg, step, t, x_faces, y_faces, U,
                _snapshot_extras(output_cfg, material_id, phi),
                materials, material_id)
            while next_output <= t + 1.0e-14:
                next_output += output_interval

    if step >= max_steps and t < final_time:
        raise RuntimeError(
            "run_case_2d: maximum step count reached before final time.")

    diagnostics.finalize(step, t, U, material_id, phi)

    return {
        "U": U,
        "phi": phi,
        "material_id": material_id,
        "x_cells": x,
        "y_cells": y,
        "x_faces": x_faces,
        "y_faces": y_faces,
        "dx": dx,
        "dy": dy,
        "t": t,
        "step": step,
        "rk_mismatch_x": last_rk_mismatch_x,
        "rk_mismatch_y": last_rk_mismatch_y,
        "parallel": parallel_info,
    }


if __name__ == "__main__":
    run()
