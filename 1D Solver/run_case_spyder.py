# Just run this file after setting the case data in the case file
# This file is the main file which connects every module in the models folder and it runs in sequence
import os
import sys
import glob
import numpy as np

# ---- make sure to run from the project root (important for Spyder case) ----
ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from case_1D_expanding_slab import get_config # reading input file data
cfg = get_config()

from models.mesh_1d import build_mesh_1d
from models.state_vars import allocate_state
from models.time_step_cfl import compute_dt_cfl
from models.rhs_euler_1d import rhs_euler_1d
from models.io_outputs import write_snapshot
from models.diagnostic_outputs_1d import DiagnosticOutput1D
from models.level_set_solution import (initialize_phi_and_material_id_1d,level_set_solution,
    find_zero_crossings_1d, compute_level_set_rhs_volume_1d,
    update_material_id_from_phi,)
from models.bi_material_cell_updates import (initialize_interfacial_flux_mismatch_1d,
    initialize_interfacial_flux_mismatch_steps_1d,
    initialize_correction_reservoir_1d,
    accumulate_interfacial_flux_mismatch_1d,
    accumulate_interfacial_flux_mismatch_steps_1d,update_bi_material_cells_1d,)
from models.bi_material_state_calculator.volume import (
    initialize_volumetric_state_var_1d, update_volumetric_state_var,
    advance_volumetric_quantities_euler_1d,
    update_volumetric_var_with_flux_1d,
    update_interface_cell_state_from_volume_1d,)

# ----------------------------
# Read run controls from input
# ----------------------------
t_final = float(cfg["time"]["max_time"])
CFL = float(cfg["time"]["cfl"])
dt_out = float(cfg["output"]["time_interval"])
out_dir = cfg["output"]["folder_name"]
sol_name = cfg["output"]["solution_name"]

# ----------------------------
# CLEAN the outputs to re-write
# ----------------------------
os.makedirs(out_dir, exist_ok=True)
pvd_path = os.path.join(out_dir, f"{sol_name}.pvd")
if os.path.exists(pvd_path):
    os.remove(pvd_path)

for f in glob.glob(os.path.join(out_dir, f"{sol_name}_*.vtk")):
    try:
        os.remove(f)
    except Exception:
        pass

for f in glob.glob(os.path.join(out_dir, f"{sol_name}_*.vtp")):
    try:
        os.remove(f)
    except Exception:
        pass

print("Writing outputs to:", os.path.abspath(out_dir))

#-----------------------------------------------------

# ----------------------------
# Building mesh # 1 time only per run
# ----------------------------
x_cells, x_faces, dx = build_mesh_1d(cfg["mesh"])
N = len(x_cells)
print("N cells:", N)

# ----------------------------
# Initial material layout / phi / U
# ----------------------------
material_defs = cfg["equations"]["materials"]
default_state = cfg["initial_condition"]["default_state"]
default_material_id = int(default_state["material_id"])
default_eos = material_defs[default_material_id]

phi, material_id = initialize_phi_and_material_id_1d(x_cells, dx, cfg) # initialize Phi

U = allocate_state(N) # lets create empty container to later fill in the values of U

vapor_region_cfg = cfg["initial_condition"].get("initial_vapor_region", {})
vapor_enabled = isinstance(vapor_region_cfg, dict) and bool(vapor_region_cfg.get("enabled", False)) #decide if there is a vapor region or not
vapor_state = vapor_region_cfg.get("state", {}) if vapor_enabled else {}
vapor_material_id = int(vapor_state.get("material_id", 1)) if vapor_enabled else 1

for i in range(N): # initilize onservative state variables U
    if vapor_enabled and material_id[i] == vapor_material_id:
        state_i = vapor_state
    else:
        state_i = default_state

    mat_i = int(state_i["material_id"])
    eos_i = material_defs[mat_i]
    nasg_i = eos_i["nasg"]

    rho_i = float(state_i["rho"])
    u_i = float(state_i["u"])
    p_i = float(state_i["p"])

    gamma_i = float(nasg_i["gamma"])
    pc_i = float(nasg_i["p_c"])
    b_i = float(nasg_i["b"])
    q_i = float(nasg_i["q"])

    V_i = 1.0 / rho_i
    e_i = q_i + (p_i + gamma_i * pc_i) * (V_i - b_i) / (gamma_i - 1.0) # from EOS
    E_i = e_i + 0.5 * u_i * u_i

    U[i, 0] = rho_i
    U[i, 1] = rho_i * u_i
    U[i, 2] = rho_i * E_i

####### Kept only for backward compatibility with some functions.
eos = default_eos



# Swept-volume method variables. They are only active when phase_change = "Volume".
volume_method = str(cfg.get("multiphase", {}).get("phase_change", "")).strip().lower() == "volume"
left_volumetric_var = None
right_volumetric_var = None
volume_interface_velocity_rk1 = None
volume_interface_velocity_rk2 = None
volume_left_flux_rk1 = None
volume_right_flux_rk1 = None
volume_left_flux_rk2 = None
volume_right_flux_rk2 = None
volume_interface_pressure_rk1 = None
volume_interface_pressure_rk2 = None
volume_moving_flux_rk1 = None
volume_moving_flux_rk2 = None
volume_left_outer_flux_rk1 = None
volume_left_outer_flux_rk2 = None
volume_right_outer_flux_rk1 = None
volume_right_outer_flux_rk2 = None
# Optional RHS diagnostics must also exist when Volume is disabled.
volume_left_star_state_rk1 = None
volume_left_star_state_rk2 = None
volume_right_star_state_rk1 = None
volume_right_star_state_rk2 = None
cell_left_flux_rk1 = None
cell_left_flux_rk2 = None
cell_right_flux_rk1 = None
cell_right_flux_rk2 = None
if volume_method:
    left_volumetric_var, right_volumetric_var = initialize_volumetric_state_var_1d(N)
    volume_interface_velocity_rk1 = np.full(N + 1, np.nan, dtype=float)
    volume_interface_velocity_rk2 = np.full(N + 1, np.nan, dtype=float)
    volume_left_flux_rk1 = np.zeros((N + 1, 3), dtype=float)
    volume_right_flux_rk1 = np.zeros((N + 1, 3), dtype=float)
    volume_left_flux_rk2 = np.zeros((N + 1, 3), dtype=float)
    volume_right_flux_rk2 = np.zeros((N + 1, 3), dtype=float)
    volume_interface_pressure_rk1 = np.full(N + 1, np.nan, dtype=float)
    volume_interface_pressure_rk2 = np.full(N + 1, np.nan, dtype=float)
    volume_moving_flux_rk1 = np.zeros((N + 1, 3), dtype=float)
    volume_moving_flux_rk2 = np.zeros((N + 1, 3), dtype=float)
    volume_left_outer_flux_rk1 = np.zeros((N, 3), dtype=float)
    volume_left_outer_flux_rk2 = np.zeros((N, 3), dtype=float)
    volume_right_outer_flux_rk1 = np.zeros((N, 3), dtype=float)
    volume_right_outer_flux_rk2 = np.zeros((N, 3), dtype=float)
    volume_left_star_state_rk1 = np.zeros((N + 1, 3), dtype=float)
    volume_left_star_state_rk2 = np.zeros((N + 1, 3), dtype=float)
    volume_right_star_state_rk1 = np.zeros((N + 1, 3), dtype=float)
    volume_right_star_state_rk2 = np.zeros((N + 1, 3), dtype=float)
    cell_left_flux_rk1 = np.zeros((N, 3), dtype=float)
    cell_right_flux_rk1 = np.zeros((N, 3), dtype=float)

    cell_left_flux_rk2 = np.zeros((N, 3), dtype=float)
    cell_right_flux_rk2 = np.zeros((N, 3), dtype=float)
    update_volumetric_state_var(
        U, phi, x_cells, dx, left_volumetric_var, right_volumetric_var)

extras = {}
if cfg["output"].get("fields", {}).get("material_id", False):
    extras["material_id"] = material_id
if cfg["output"].get("fields", {}).get("levelset0", False):
    extras["levelset0"] = phi

# ----------------------------
# Numerics | Need to update it later
# ----------------------------
numerics_cfg = cfg['space']['navier_stokes'] # Flux and reconstruction choices from the case file
bc_cfg = cfg

# mismatch and two temporary RK-stage arrays.
accumulated_flux_mismatch = initialize_interfacial_flux_mismatch_1d(N)
accumulated_flux_mismatch_steps = initialize_interfacial_flux_mismatch_steps_1d(N)
correction_reservoir, correction_steps_remaining = (
    initialize_correction_reservoir_1d(N))
flux_mismatch_rk1 = np.zeros_like(accumulated_flux_mismatch)
flux_mismatch_rk2 = np.zeros_like(accumulated_flux_mismatch)

# ----------------------------
# Time loop + output
# ----------------------------
t = 0.0
step = 0
max_steps = 100000 # safety limit; normal run stops at t_final
next_out = 0.0

write_snapshot(
    cfg["output"], step, t, x_cells, U, eos=eos, extras=extras,
    materials=material_defs, material_id=material_id,)


diagnostics = DiagnosticOutput1D(cfg["output"], x_cells, dx, material_defs)
if volume_method:
    diagnostics.initialize(
        step, t, U, material_id, phi,
        left_volumetric_var, right_volumetric_var)
else:
    diagnostics.initialize(step, t, U, material_id, phi)

next_out += dt_out

print("Start time loop: t_final=", t_final, " CFL=", CFL, " dt_out=", dt_out)

while t < t_final and step < max_steps:
    dt = compute_dt_cfl(U, dx, material_defs, material_id, CFL)
    if t + dt > t_final: # for the last time step
        dt = t_final - t

    # RK2 update for Euler state
    U_old = U.copy()
    
    #diagnosis
    if volume_method:
        interface_positions_old = find_zero_crossings_1d(x_cells, phi)
        x_interface_old = float(interface_positions_old[0])

        diagnostic_i = int(
            np.searchsorted(x_faces, x_interface_old, side="right") - 1)

        diagnostic_cells = [
            diagnostic_i - 1,
            diagnostic_i,
            diagnostic_i + 1,
        ]

        print("\n============================================================")
        print(f"DIAGNOSTIC TIME STEP {step + 1}")
        print(f"t_n = {t:.12e}, dt = {dt:.12e}")
        print(
            f"interface at t_n = {x_interface_old:.12e}, "
            f"interface cell i = {diagnostic_i}")

        print("INITIAL STATES U at t_n [rho, rho*u, rho*E]")

        for j in diagnostic_cells:
            print(
                f"cell {j}: x = {x_cells[j]:.12e}, "
                f"material_id = {material_id[j]}, "
                f"phi = {phi[j]:.12e}, "
                f"U = {U_old[j, :]}")

    #diagnosis ends
    
    interface_callback_rk1 = diagnostics.make_interface_callback(step + 1, t, "RK1")
    k1 = rhs_euler_1d(
        U, x_faces, dx, material_defs, material_id, numerics_cfg, bc_cfg,
        interface_callback=interface_callback_rk1,
        flux_mismatch=flux_mismatch_rk1,
        x_cells=x_cells,
        phi=phi,
        left_volumetric_var=left_volumetric_var,
        right_volumetric_var=right_volumetric_var,
        volume_interface_velocity=volume_interface_velocity_rk1,
        volume_interface_pressure=volume_interface_pressure_rk1,
        volume_moving_flux_output=volume_moving_flux_rk1,
        volume_left_outer_flux_output=volume_left_outer_flux_rk1,
        volume_right_outer_flux_output=volume_right_outer_flux_rk1,
        volume_left_star_state_output=volume_left_star_state_rk1,
        volume_right_star_state_output=volume_right_star_state_rk1,
        volume_left_flux_output=volume_left_flux_rk1,
        volume_right_flux_output=volume_right_flux_rk1,
        cell_left_flux_output=cell_left_flux_rk1,
        cell_right_flux_output=cell_right_flux_rk1)
    
    U1_euler = U + dt * k1

    if volume_method:
        # RK1 predictor for the moving interface and its two sub-cell Q values.
        phi_stage = phi + dt * compute_level_set_rhs_volume_1d(
            U, phi, x_cells, dx, volume_interface_velocity_rk1)

        left_volumetric_stage, right_volumetric_stage = (
            advance_volumetric_quantities_euler_1d(
                U, U1_euler, phi, phi_stage, x_cells, dx, dt,
                left_volumetric_var, right_volumetric_var,
                volume_interface_velocity_rk1, volume_moving_flux_rk1,
                volume_left_outer_flux_rk1, volume_right_outer_flux_rk1,
                volume_left_flux_rk1, volume_right_flux_rk1,
                volume_left_star_state_rk1, volume_right_star_state_rk1,
                cell_left_flux_rk1, cell_right_flux_rk1,
                material_defs, material_id))

        # The RK2 Euler stage must use the state represented by the RK1
        # volumetric quantities and the RK1 interface position.
        U1 = update_interface_cell_state_from_volume_1d(
            U1_euler, phi_stage, x_cells, dx,
            left_volumetric_stage, right_volumetric_stage)

        material_id_stage = update_material_id_from_phi(
            phi_stage, material_id, default_material_id, vapor_material_id)
    else:
        U1 = U1_euler
        phi_stage = phi
        material_id_stage = material_id
        left_volumetric_stage = left_volumetric_var
        right_volumetric_stage = right_volumetric_var
    accumulate_interfacial_flux_mismatch_1d(accumulated_flux_mismatch, flux_mismatch_rk1, 0.5 * dt)

    interface_callback_rk2 = diagnostics.make_interface_callback(step + 1, t + dt, "RK2")
    k2 = rhs_euler_1d(
        U1, x_faces, dx, material_defs, material_id_stage, numerics_cfg, bc_cfg,
        interface_callback=interface_callback_rk2,
        flux_mismatch=flux_mismatch_rk2,
        x_cells=x_cells,
        phi=phi_stage,
        left_volumetric_var=left_volumetric_stage,
        right_volumetric_var=right_volumetric_stage,
        volume_interface_velocity=volume_interface_velocity_rk2,
        volume_interface_pressure=volume_interface_pressure_rk2,
        volume_moving_flux_output=volume_moving_flux_rk2,
        volume_left_outer_flux_output=volume_left_outer_flux_rk2,
        volume_right_outer_flux_output=volume_right_outer_flux_rk2,
        volume_left_star_state_output=volume_left_star_state_rk2,
        volume_right_star_state_output=volume_right_star_state_rk2,
        volume_left_flux_output=volume_left_flux_rk2,
        volume_right_flux_output=volume_right_flux_rk2,
        cell_left_flux_output=cell_left_flux_rk2,
        cell_right_flux_output=cell_right_flux_rk2)
    
    
    U_candidate = 0.5 * U + 0.5 * (U1_euler + dt * k2)
    accumulate_interfacial_flux_mismatch_1d(accumulated_flux_mismatch, flux_mismatch_rk2, 0.5 * dt)
    accumulate_interfacial_flux_mismatch_steps_1d(
        accumulated_flux_mismatch_steps, flux_mismatch_rk1, flux_mismatch_rk2)
    
    if volume_method:
        active_faces_rk1 = np.where(
            np.isfinite(volume_interface_velocity_rk1))[0]
        active_faces_rk2 = np.where(
            np.isfinite(volume_interface_velocity_rk2))[0]

        interface_face_rk1 = int(active_faces_rk1[0])
        interface_face_rk2 = int(active_faces_rk2[0])

        print("RK1 AND RK2 FLUXES USED BY CELLS")

        for j in diagnostic_cells:
            print(
                f"cell {j}: "
                f"RK1 left = {cell_left_flux_rk1[j, :]}, "
                f"RK1 right = {cell_right_flux_rk1[j, :]}")

            print(
                f"         "
                f"RK2 left = {cell_left_flux_rk2[j, :]}, "
                f"RK2 right = {cell_right_flux_rk2[j, :]}")

        print(
            f"bi-material face RK1 = {interface_face_rk1}, "
            f"x_face = {x_faces[interface_face_rk1]:.12e}")
        print(
            f"bi-material face RK2 = {interface_face_rk2}, "
            f"x_face = {x_faces[interface_face_rk2]:.12e}")

        print(
            f"bi-material FL RK1 = "
            f"{volume_left_flux_rk1[interface_face_rk1, :]}")

        print(
            f"bi-material FR RK1 = "
            f"{volume_right_flux_rk1[interface_face_rk1, :]}")

        print(
            f"bi-material FL RK2 = "
            f"{volume_left_flux_rk2[interface_face_rk2, :]}")

        print(
            f"bi-material FR RK2 = "
            f"{volume_right_flux_rk2[interface_face_rk2, :]}")

        print(
            f"V* RK1 = "
            f"{volume_interface_velocity_rk1[interface_face_rk1]:.12e}")

        print(
            f"V* RK2 = "
            f"{volume_interface_velocity_rk2[interface_face_rk2]:.12e}")

        active_pressure_rk1 = volume_interface_pressure_rk1[
            np.isfinite(volume_interface_pressure_rk1)][0]
        active_pressure_rk2 = volume_interface_pressure_rk2[
            np.isfinite(volume_interface_pressure_rk2)][0]
        active_moving_flux_rk1 = volume_moving_flux_rk1[
            np.isfinite(volume_interface_pressure_rk1), :][0]
        active_moving_flux_rk2 = volume_moving_flux_rk2[
            np.isfinite(volume_interface_pressure_rk2), :][0]

        print(f"p* RK1 = {active_pressure_rk1:.12e}")
        print(f"p* RK2 = {active_pressure_rk2:.12e}")
        print(f"moving interface G RK1 = {active_moving_flux_rk1}")
        print(f"moving interface G RK2 = {active_moving_flux_rk2}")

    # Level-set update and material-ID update
    material_id_old = material_id.copy()
    phi_old = phi.copy()
    phi, material_id_new, _ = level_set_solution(U,U1,phi,material_id,x_cells,dx,dt,step + 1,cfg,
        volume_interface_velocity_rk1, volume_interface_velocity_rk2,)

    if volume_method:
        
        interface_positions_new = find_zero_crossings_1d(x_cells, phi)
        x_interface_new = float(interface_positions_new[0])

        diagnostic_i_new = int(
            np.searchsorted(x_faces, x_interface_new, side="right") - 1)

        left_volume_old = (
            x_interface_old - x_faces[diagnostic_i])

        right_volume_old = (
            x_faces[diagnostic_i + 1] - x_interface_old)

        left_volume_new = (
            x_interface_new - x_faces[diagnostic_i_new])

        right_volume_new = (
            x_faces[diagnostic_i_new + 1] - x_interface_new)

        left_accumulated_flux_quantity = 0.5 * dt * (
            volume_left_flux_rk1[interface_face_rk1, :]
            + volume_left_flux_rk2[interface_face_rk2, :])

        right_accumulated_flux_quantity = 0.5 * dt * (
            volume_right_flux_rk1[interface_face_rk1, :]
            + volume_right_flux_rk2[interface_face_rk2, :])

        print("INTERFACE MOVEMENT / VOLUMES")

        print(f"interface at t_n     = {x_interface_old:.12e}")
        print(f"interface at t_n+1   = {x_interface_new:.12e}")

        print(
            f"delta x = delta V    = "
            f"{x_interface_new - x_interface_old:.12e}")

        print(f"old interface cell   = {diagnostic_i}")
        print(f"new interface cell   = {diagnostic_i_new}")

        print(f"left volume at t_n   = {left_volume_old:.12e}")
        print(f"right volume at t_n  = {right_volume_old:.12e}")

        print(f"left volume at t_n+1 = {left_volume_new:.12e}")
        print(f"right volume at t_n+1= {right_volume_new:.12e}")

        print(
            f"left fixed-grid bi-material flux quantity  = "
            f"{left_accumulated_flux_quantity}")

        print(
            f"right fixed-grid bi-material flux quantity = "
            f"{right_accumulated_flux_quantity}")

        moving_interface_quantity = 0.5 * dt * (
            active_moving_flux_rk1 + active_moving_flux_rk2)

        print(
            f"moving-interface quantity G*dt = "
            f"{moving_interface_quantity}")
        
        
        # Update QL and QR with the moving-control-volume equations:
        # dQL/dt = F_left_outer - G,
        # dQR/dt = G - F_right_outer,
        # where G = [0, p*, p*V*].
        volume_delta, volume_direction = update_volumetric_var_with_flux_1d(
            U_old, U_candidate, phi_old, phi_stage, phi,
            x_cells, dx, dt,
            left_volumetric_var, right_volumetric_var,
            left_volumetric_stage, right_volumetric_stage,
            volume_interface_velocity_rk1, volume_interface_velocity_rk2,
            volume_moving_flux_rk1, volume_moving_flux_rk2,
            volume_left_outer_flux_rk1, volume_left_outer_flux_rk2,
            volume_right_outer_flux_rk1, volume_right_outer_flux_rk2,
            volume_left_flux_rk1, volume_left_flux_rk2,
            volume_right_flux_rk1, volume_right_flux_rk2,
            volume_left_star_state_rk1, volume_left_star_state_rk2,
            volume_right_star_state_rk1, volume_right_star_state_rk2,
            cell_left_flux_rk1, cell_left_flux_rk2,
            cell_right_flux_rk1, cell_right_flux_rk2,
            material_defs, material_id_old)

        # Old Q_swept neighbor transfer is intentionally not used.
        volume_swept_quantity = np.zeros((1, 3), dtype=float)
        
        
        print("RUN_CASE: VOLUMETRIC VARIABLES AFTER COMPLETE VOLUME UPDATE")
        print(f"time step = {step + 1}")
        print(f"interface cell i = {diagnostic_i_new}")
        print(
    f"left_volumetric_var[{diagnostic_i_new}]  = "
    f"{left_volumetric_var[diagnostic_i_new, :]}")
        print(
    f"right_volumetric_var[{diagnostic_i_new}] = "
    f"{right_volumetric_var[diagnostic_i_new, :]}")

    # Update conservative states after the new material IDs are known.
    U_new, crossed_cells = update_bi_material_cells_1d(U_candidate, material_id_old,
    material_id_new, x_cells, dx, cfg, accumulated_flux_mismatch,
    accumulated_flux_mismatch_steps, correction_reservoir,
    correction_steps_remaining, phi, left_volumetric_var,
    right_volumetric_var)

    material_id = material_id_new
    U = U_new
    
    
    if volume_method:
        print("UPDATED STATES U at t_n+1 [rho, rho*u, rho*E]")

        for j in diagnostic_cells:
            print(
                f"cell {j}: x = {x_cells[j]:.12e}, "
                f"material_id = {material_id[j]}, "
                f"phi = {phi[j]:.12e}, "
                f"U = {U[j, :]}")

        print(
            f"left_volumetric_var[{diagnostic_i_new}]  = "
            f"{left_volumetric_var[diagnostic_i_new, :]}")

        print(
            f"right_volumetric_var[{diagnostic_i_new}] = "
            f"{right_volumetric_var[diagnostic_i_new, :]}")

        print("============================================================")


    if cfg["output"].get("fields", {}).get("material_id", False):
        extras["material_id"] = material_id
    if cfg["output"].get("fields", {}).get("levelset0", False):
        extras["levelset0"] = phi
    if cfg["output"].get("debug", False):
        extras["crossed_cells"] = crossed_cells.astype(float)


    t += dt
    step += 1
    print(f"time: step={step}, t={t:.6e},")
    
    #interface_positions = find_zero_crossings_1d(x_cells, phi)

    #print(
        #f"Final interface position: step={step}, t={t:.6e}, "
        #f"phi=0 at x={interface_positions}")
    
    
    
    if volume_method:
        diagnostics.record_step(
            step, t, dt, U_old, U_candidate, U_new, material_id_old,
            material_id, crossed_cells, phi,
            left_volumetric_var, right_volumetric_var)
    else:
        diagnostics.record_step(
            step, t, dt, U_old, U_candidate, U_new, material_id_old,
            material_id, crossed_cells, phi)

    if t >= next_out - 1e-15:
        write_snapshot(cfg["output"],step, t, x_cells, U,eos=eos,extras=extras,
            materials=material_defs,material_id=material_id,)
        
        next_out += dt_out
        print(f"WROTE: step={step}, t={t:.6e}")

if volume_method:
    diagnostics.finalize(
        step, t, U, material_id, phi,
        left_volumetric_var, right_volumetric_var)
else:
    diagnostics.finalize(step, t, U, material_id, phi)

print("DONE. final step =", step, "final t =", t)
print("Open in ParaView:", os.path.join(os.path.abspath(out_dir), f"{sol_name}.pvd"))