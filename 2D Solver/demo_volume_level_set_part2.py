"""Standalone PART-2 demonstration: prescribed Riemann star velocities.

Run: python demo_volume_level_set_part2.py

This tests level-set geometry/velocity-extension/RK2 only.  Full Volume
sub-volume updates and Euler coupling are reserved for Parts 3/4.
"""

import copy
import numpy as np

from case_2D_CircleAdvection import get_config
from models.mesh_2d import build_mesh_2d
from models.level_set_solution import initialize_phi_and_material_id_2d, level_set_solution
from models.level_set_volume_2d import compute_level_set_rhs_volume_2d
from models.interface_reconstruction import reconstruct_interface_2d


def velocity_provider(_phi, segments, _stage):
    """Demo only: identical uniform material velocities imply this star velocity."""
    return np.tile([0.45, 0.25], (len(segments), 1))


def main():
    cfg = copy.deepcopy(get_config())
    cfg["mesh"].update({"nx": 80, "ny": 80})
    cfg["multiphase"]["phase_change"] = "Volume"
    cfg["multiphase"]["volume_level_set"] = {
        "core_cells": 5,
        "transition_cells": 10,
        "node_band_radius": 5,
    }
    cfg["space"]["level_set"][0]["reinitialization"]["frequency"] = 0

    x, x_faces, dx, y, y_faces, dy = build_mesh_2d(cfg["mesh"])
    phi, material_id = initialize_phi_and_material_id_2d(x, y, dx, dy, cfg)
    phi_initial = phi.copy()
    U = np.zeros(phi.shape + (4,))
    U[..., 0] = 1.0
    U[..., 1] = 0.45
    U[..., 2] = 0.25
    U[..., 3] = 3.0

    initial_geometry = reconstruct_interface_2d(phi, x_faces, y_faces)
    rhs, fields = compute_level_set_rhs_volume_2d(
        U, phi, x, y, dx, dy, velocity_provider, cfg, return_fields=True)

    dt = 0.002
    nsteps = 24
    for step in range(1, nsteps + 1):
        phi, material_id, _ = level_set_solution(
            U, U, phi, material_id, x, y, dx, dy,
            dt, step, cfg, volume_velocity_provider=velocity_provider)
    end_geometry = reconstruct_interface_2d(phi, x_faces, y_faces)

    X, Y = np.meshgrid(x, y)
    area = end_geometry["area_1"]
    centroid_x = float(np.sum(area * X) / np.sum(area))
    centroid_y = float(np.sum(area * Y) / np.sum(area))
    expected_x = 0.35 + 0.45 * dt * nsteps
    expected_y = 0.40 + 0.25 * dt * nsteps
    print(f"Prescribed uniform velocity: (0.45, 0.25)")
    print(f"Initial interface segments: {len(initial_geometry['segments'])}")
    print(f"Expected circle center: ({expected_x:.6f}, {expected_y:.6f})")
    print(f"Measured area-weighted center: ({centroid_x:.6f}, {centroid_y:.6f})")
    print(f"Center errors: ({centroid_x-expected_x:+.3e}, {centroid_y-expected_y:+.3e})")
    print(f"Initial area: {np.sum(initial_geometry['area_1']):.8f}")
    print(f"Final area:   {np.sum(end_geometry['area_1']):.8f}")

    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        axes[0].contour(X, Y, phi_initial, [0], linewidths=2)
        axes[0].contour(X, Y, phi, [0], linestyles="--", linewidths=2)
        axes[0].plot(expected_x, expected_y, marker="x", markersize=9)
        axes[0].set_title("Initial and translated zero contour")
        axes[0].set_xlabel("x")
        axes[0].set_ylabel("y")
        axes[0].set_aspect("equal")
        axes[1].imshow(fields["weight"], origin="lower", extent=[0, 1, 0, 1],
                       vmin=0, vmax=1)
        axes[1].contour(X, Y, phi_initial, [0])
        axes[1].set_title("Core + transition weight")
        axes[1].set_xlabel("x")
        axes[1].set_ylabel("y")
        axes[1].set_aspect("equal")
        fig.savefig("volume_level_set_part2_demo.png", dpi=150)
        plt.close(fig)
        print("Saved: volume_level_set_part2_demo.png")
    except ImportError:
        print("Matplotlib not found: numerical demo succeeded, plot omitted.")


if __name__ == "__main__":
    main()
