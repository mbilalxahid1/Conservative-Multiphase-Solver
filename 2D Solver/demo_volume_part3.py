"""Independent Part-3 demonstration of the 2D conservative Volume method.

Run from the project directory:
    python demo_volume_part3.py

Uses the real interface reconstruction, two-material Riemann solver,
normal-speed extension, aperture-weighted face fluxes, SSP-RK2 and
conservative local crossing/merging. The *main* run_case_2d.py remains on
its original FIVER/Extrapolation path until the next integration step.
"""

import copy

import numpy as np

from case_2D_CircleAdvection import get_config
from models.mesh_2d import build_mesh_2d
from models.geometry_2d import initial_state_at_point_2d
from models.level_set_solution import initialize_phi_and_material_id_2d
from models.eos_nasg import primitive_to_conservative_2d
from models.bi_material_state_calculator.volume import (
    initialize_volume_state_2d, reconstruct_volume_geometry_2d,
    volume_step_2d,
)


def main():
    cfg = copy.deepcopy(get_config())
    cfg["mesh"]["nx"] = 32
    cfg["mesh"]["ny"] = 32
    cfg["parallel"]["enabled"] = False
    cfg["multiphase"]["phase_change"] = "Volume"
    cfg["multiphase"]["volume_solver"] = {
        "merge_area_fraction": 0.001,
        "max_hops": 2,
        "max_interface_cfl": 0.5,
    }
    cfg["multiphase"]["volume_level_set"] = {
        "core_cells": 5,
        "transition_cells": 20,
        "node_band_radius": 5,
    }
    cfg["space"]["level_set"][0]["reinitialization"]["frequency"] = 0
    for state in (cfg["initial_condition"]["default_state"],
                  cfg["initial_condition"]["geometry"]["state"]):
        state.update({"u": 0.1, "v": 0.05, "p": 1.0})

    x, xf, dx, y, yf, dy = build_mesh_2d(cfg["mesh"])
    phi, ids = initialize_phi_and_material_id_2d(x, y, dx, dy, cfg)
    phi_start = phi.copy()
    materials = cfg["equations"]["materials"]
    U = np.zeros(phi.shape + (4,))
    for j in range(len(y)):
        for i in range(len(x)):
            state = initial_state_at_point_2d(x[i], y[j], cfg)
            U[j, i] = primitive_to_conservative_2d(
                state["rho"], state["u"], state["v"],
                state["p"], materials[state["material_id"]])
    geometry = reconstruct_volume_geometry_2d(phi, xf, yf, cfg)
    Q0, Q1 = initialize_volume_state_2d(U, ids, geometry, materials, x, y)

    mass0_initial = Q0[..., 0].sum()
    mass1_initial = Q1[..., 0].sum()
    area_initial = geometry["area_1"].sum()
    dt = 0.003
    steps = 8
    conservation_errors = []
    for step in range(1, steps + 1):
        Q0, Q1, phi, ids, info = volume_step_2d(
            Q0, Q1, phi, xf, yf, cfg, dt, step=step)
        conservation_errors.append(np.max(abs(info["conservation_error"])))
    final_geometry = reconstruct_volume_geometry_2d(phi, xf, yf, cfg)
    XX, YY = np.meshgrid(x, y)
    A = final_geometry["area_1"]
    measured_x = np.sum(A*XX)/np.sum(A)
    measured_y = np.sum(A*YY)/np.sum(A)
    initial_center = cfg["initial_condition"]["geometry"]
    target_x = initial_center["center_x"] + 0.1*dt*steps
    target_y = initial_center["center_y"] + 0.05*dt*steps
    print("2D Volume independent RK2 demonstration")
    print("  initial material-0 mass:", mass0_initial)
    print("  final material-0 mass:  ", Q0[..., 0].sum())
    print("  initial material-1 mass:", mass1_initial)
    print("  final material-1 mass:  ", Q1[..., 0].sum())
    print("  material-1 area change: ", A.sum()-area_initial)
    print("  exact translation center:", (target_x, target_y))
    print("  observed area centroid:  ", (measured_x, measured_y))
    print("  max stage budget error: ", max(conservation_errors))

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    axes[0].contour(XX, YY, phi_start, [0.0], colors="tab:blue", linewidths=2)
    axes[0].contour(XX, YY, phi, [0.0], colors="tab:orange",
                    linestyles="--", linewidths=2)
    axes[0].scatter([target_x], [target_y], marker="x", color="black", s=50)
    axes[0].set_title("Initial (blue), final (dashed), target center")
    axes[1].imshow(final_geometry["area_fraction_1"], origin="lower",
                   extent=(xf[0], xf[-1], yf[0], yf[-1]),
                   vmin=0, vmax=1, cmap="Blues")
    axes[1].set_title("Final reconstructed material-1 area fraction")
    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
    fig.savefig("volume_part3_circle_demo.png", dpi=160)
    plt.close(fig)
    print("  output: volume_part3_circle_demo.png")


if __name__ == "__main__":
    main()
