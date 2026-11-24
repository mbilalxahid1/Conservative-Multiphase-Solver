# 2D shock-induced collapse of an air bubble in water.
#
# Benchmark sources:
#   Maltsev, Skote & Tsoutsanis, JCP 502 (2024) 112819, Sec. 3.5.
#   Fechter et al., Journal of Scientific Computing (2023), Sec. 5.3.
#
# Published geometry/states:
#   domain:       [0, 0.024] x [0, 0.024] m
#   bubble:       center (0.012, 0.012) m, radius 0.003 m
#   shock front:  x = 0.0066 m, moving right
#   shock Mach:   1.72
#
# IMPORTANT:
# The circular bubble is the only material interface and is tracked by phi.
# The planar shock is represented by initial_condition.state_regions and is
# only a same-material primitive-state discontinuity in water.

CONFIG = {
    "mesh": {
        "type": "TwoDimensional",
        "x0": 0.0,
        "xmax": 0.024,
        "y0": 0.0,
        "ymax": 0.024,

        # D_b = 0.006 m.  320 x 320 gives 80 cells across the initial bubble
        # diameter, matching the approximate resolution quoted by Maltsev et al.
        "use_uniform_x": True,
        "nx": 520,
        "use_control_points_x": False,
        "control_points_x": [],

        "use_uniform_y": True,
        "ny": 520,
        "use_control_points_y": False,
        "control_points_y": [],

        # Published: non-reflecting left/right, Euler-slip top/bottom.
        # Current Python 'outlet' is a zero-gradient approximation to a
        # non-reflecting boundary; wall is the Euler-slip wall.
        "bc_x0": "outlet",
        "bc_xmax": "outlet",
        "bc_y0": "wall",
        "bc_ymax": "wall",
    },

    "boundary_conditions": {
        "outlet": {"type": "outlet"},
        "wall": {"type": "wall"},
    },

    "equations": {
        "materials": {
            0: {
                "name": "water",
                "eos": "NobleAbelStiffenedGas",
                # Stiffened gas = NASG with b=q=0.
                "nasg": {
                    "gamma": 4.4,
                    "p_c": 6.0e8,
                    "b": 0.0,
                    "q": 0.0,
                    "cv": 1.0,
                },
            },
            1: {
                "name": "air_bubble",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4,
                    "p_c": 0.0,
                    "b": 0.0,
                    "q": 0.0,
                    "cv": 1.0,
                },
            },
        },
        "material_transitions": [],
    },

    "initial_condition": {
        # Pre-shock water fills the background.
        "default_state": {
            "material_id": 0,
            "rho": 1000.0,
            "u": 0.0,
            "v": 0.0,
            "p": 1.0e5,
        },

        # Same water material, different primitive state.  This is NOT a
        # material interface and therefore does not create another level set.
        "state_regions": [{
            "type": "Rectangle",
            "x_min": 0.0,
            "x_max": 0.0066,
            "y_min": 0.0,
            "y_max": 0.024,
            "state": {
                "material_id": 0,
                "rho": 1323.65,
                "u": 681.058,
                "v": 0.0,
                "p": 1.9e9,
            },
        }],

        # Air/water material interface.
        "geometry": {
            "enabled": True,
            "type": "Circle",
            "center_x": 0.012,
            "center_y": 0.012,
            "radius": 0.003,
            "state": {
                "material_id": 1,
                "rho": 1.0,
                "u": 0.0,
                "v": 0.0,
                "p": 1.0e5,
            },
        },
    },

    "space": {
        "navier_stokes": {
            "flux": "Rusanov",
            "reconstruction": {"type": "Constant"},
        },
        "level_set": [{
            "material_id": 1,
            # Reinitialize after each accepted RK2 step, matching the
            # shock/bubble settings supplied for this calculation.
            "reinitialization": {
                "frequency": 1,
                "max_iterations": 40,
                "cfl": 0.4,
                "convergence_tolerance": 1.0e-7,
            },
            "reconstruction": {
                # Select the level-set spatial derivative (both SSP-RK2 stages).
                # Options: "FirstOrderUpwind", "SecondOrderENO",
                # or "ThirdOrderUpwindFD" (original scheme).
                "type": "ThirdOrderUpwindFD",
                "limiter": "None",
            },
        }],
    },

    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        "phase_change": "Volume",
        "phase_change_direction": "Upwind",
        "volume_solver": {
            # "Directional": original independent X/Y + selected-normal solve.
            # "SingleInterfaceNormal": ALWAYS ONE geometric-normal Riemann
            # solve for interface momentum, pressure work and level-set motion.
            "interface_riemann_method": "SingleInterfaceNormal",
            "normal_mode": "Interface",  # legacy setting for Directional only
            # Used ONLY by the original "Directional" method, if an X/Y solve
            # fails. SingleInterfaceNormal already uses a single solve always.
            "normal_star_flux_fallback": False,
            "riemann_fallback_debug": True,
            # All other Volume options retain their original defaults.
            "failure_diagnostics": True,
            "failure_debug_radius": 1,
            # First use the legacy two-cell merge and halve dt three times.
            # On the fourth attempt, enable a conservative connected-cell
            # merge of up to 5 cells if the ordinary two-cell merge fails.
            # Remaining unsuccessful attempts continue halving dt.
            "max_retries": 10,
            "adaptive_merge_enabled": True,
            "adaptive_merge_after_retries": 3,
            "adaptive_merge_max_cells": 5,
            # Replace the historical post-RK 1% Riemann-star stabilization
            # with an on-failure PRE-Riemann conservative connected-cell pool.
            # Retry dt0, dt0/2, dt0/4 without this fallback; on retry 3
            # activate it at dt0/8. Continue halving through retry 10.
            # NO material area-fraction threshold for pre-Riemann recovery.
            "small_cell_fraction": 0.0,
            "pre_riemann_merge_enabled": True,
            "pre_riemann_merge_after_retries": 3,
            "pre_riemann_merge_max_cells": 5,
            "pre_riemann_merge_max_events_per_stage": 8,
            # Fully disappeared orphan sub-volumes with no connected same-
            # material recipient are cleared automatically in the Volume
            # solver (M2C-inspired). There is deliberately NO input switch.
            # Accepted losses are logged in volume_orphan_discards.txt.
            # Only after three failed trials (dt0, dt0/2, dt0/4),
            # first try the two-cell/adaptive merge, then extrapolate an
            # EOS-invalid state from valid face-connected SAME-material cells.
            # Its Q correction is compensated in the donors and limited to
            # keep every affected cell admissible; otherwise halve dt again.
            "conservative_extrapolation_enabled": True,
            "conservative_extrapolation_after_retries": 2,
            "conservative_extrapolation_max_hops": 3,
            "conservative_extrapolation_max_donors": 8,
            "conservative_extrapolation_pressure_margin_ratio": 1.e-6,
        },
    },

    "parallel": {
        "enabled": True,
        "backend": "numba",
        "threads": 0,
        "min_cells": 2000,
    },

    "time": {
        "type": "Explicit",
        "max_time": 6.2e-6,
        "cfl": 0.4,
        "explicit": {"type": "RungeKutta2"},
    },

    "output": {
        "folder_name": "results_shock_bubble_collapse_2d/",
        "prefix": "results_shock_bubble_collapse_2d/",
        "solution_name": "solution",

        # 0.03 us output interval. Published bubble shape snapshots are at
        # 3.0, 4.1, 4.7 and 6.0 us.
        "time_interval": 2e-8,

        "fields": {
            "density": True,
            "velocity": True,
            "pressure": True,
            "temperature": False,
            "sound_speed": False,
            "material_id": True,
            "levelset0": True,
        },

        # Maltsev et al. compare centerline pressure at about 5 us.
        "line_plots": [{
            "enabled": True,
            "file_name": "centerline_y_0p012",
            "num_points": 640,
            "x0": 0.0,
            "xmax": 0.024,
            "y0": 0.012,
            "ymax": 0.012,
            "time_interval": 1e-7,
        }],

        "integrations": [{
            "enabled": True,
            "every_time_step": False,
            "time_interval": 1e-8,
            "volume_file": "volume.txt",
            "mass_file": "mass.txt",
            "momentum_file": "momentum.txt",
            "total_energy_file": "totalenergy.txt",
            "internal_energy_file": "internalenergy.txt",
        }],

        "log": {
            "enabled": True,
            "file_name": "run.log",
            "log_every_step": False,
            "material_changes": False,
            "material_change_states": False,
            "interface_fluxes": False,
            "interface_flux_stages": ["RK1", "RK2"],
        },
    },
}


def get_config():
    return CONFIG
