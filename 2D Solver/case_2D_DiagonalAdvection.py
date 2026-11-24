# Genuine 2D material-advection test.
# A square material-1 patch is transported diagonally with u>0 and v>0.
# Pressure and velocity are initially continuous across the interface, so the
# exact solution is a translating material contact. This case exercises both
# x-normal and y-normal FIVER faces, level-set motion in x and y, crossed-cell
# extrapolation, and periodic reinitialization without introducing a shock.

CONFIG = {
    "mesh": {
        "type": "TwoDimensional",
        "x0": 0.0, "xmax": 1.0,
        "y0": 0.0, "ymax": 1.0,
        "use_uniform_x": True,
        "nx": 20,
        "use_control_points_x": False,
        "control_points_x": [],
        "use_uniform_y": True,
        "ny": 20,
        "use_control_points_y": False,
        "control_points_y": [],
        "bc_x0": "outlet",
        "bc_xmax": "outlet",
        "bc_y0": "outlet",
        "bc_ymax": "outlet",
    },

    "boundary_conditions": {
        "outlet": {"type": "outlet"},
    },

    "equations": {
        "materials": {
            0: {
                "name": "background",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4, "p_c": 0.0,
                    "b": 0.0, "q": 0.0, "cv": 1.0,
                },
            },
            1: {
                "name": "patch",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4, "p_c": 0.0,
                    "b": 0.0, "q": 0.0, "cv": 1.0,
                },
            },
        },
        "material_transitions": [],
    },

    "initial_condition": {
        "default_state": {
            "material_id": 0,
            "rho": 1.0,
            "u": 40.0,
            "v": 30.0,
            "p": 1.0e5,
        },
        "geometry": {
            "enabled": True,
            "type": "Square",
            "center_x": 0.425,
            "center_y": 0.425,
            "side_length": 0.25,
            "state": {
                "material_id": 1,
                "rho": 2.0,
                "u": 40.0,
                "v": 30.0,
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
            "reinitialization": {
                "frequency": 5,
                "max_iterations": 30,
                "cfl": 0.6,
                "convergence_tolerance": 1.0e-3,
            },
            "reconstruction": {
                "type": "ThirdOrderUpwindFD",
                "limiter": "None",
            },
        }],
    },

    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        "phase_change": "Extrapolation",
        "phase_change_direction": "Upwind",
    },

    # -------------------------------------------------------------------------
    # Shared-memory CPU parallelization
    # -------------------------------------------------------------------------
    "parallel": {
        "enabled": True,
        "backend": "numba",
        # 0 = use Numba's default thread count (normally all available cores).
        # Set e.g. 4 or 8 for a fixed number of CPU threads.
        "threads": 0,
        # Small grids use the readable serial path because thread overhead can
        # dominate. Larger cases automatically use the compiled Numba kernels.
        "min_cells": 2000,
    },

    "time": {
        "type": "Explicit",
        "max_time": 1.5e-3,
        "cfl": 0.35,
        "explicit": {"type": "RungeKutta2"},
    },

    "output": {
        "folder_name": "results_diagonal_advection_2d/",
        "prefix": "results_diagonal_advection_2d/",
        "solution_name": "solution",
        "time_interval": 3.0e-4,
        "fields": {
            "density": True,
            "velocity": True,
            "pressure": True,
            "temperature": False,
            "sound_speed": False,
            "material_id": True,
            "levelset0": True,
        },
        "line_plots": [{
            "enabled": True,
            "file_name": "line_diagonal",
            "num_points": 80,
            "x0": 0.0, "xmax": 1.0,
            "y0": 0.0, "ymax": 1.0,
            "time_interval": 3.0e-4,
        }],
        "integrations": [{
            "enabled": True,
            "every_time_step": True,
            "volume_file": "volume.txt",
            "mass_file": "mass.txt",
            "momentum_file": "momentum.txt",
            "total_energy_file": "totalenergy.txt",
            "internal_energy_file": "internalenergy.txt",
        }],
        "log": {
            "enabled": True,
            "file_name": "run.log",
            "log_every_step": True,
            "material_changes": True,
            "material_change_states": True,
            "interface_fluxes": False,
            "interface_flux_stages": ["RK1", "RK2"],
        },
    },
}


def get_config():
    return CONFIG
