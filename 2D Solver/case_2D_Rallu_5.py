# 2D extrusion of Arthur Rallu thesis Fig. 4.15 shock-tube test.
# The physical solution is initially uniform in y; ny=5 is used to test the 2D code path.
# Units: kg, m, s, Pa.

CONFIG = {
    "mesh": {
        "type": "TwoDimensional",
        "x0": 0.0, "xmax": 1.0,
        "y0": 0.0, "ymax": 0.05,

        # x direction: uniform example used in this case
        "use_uniform_x": True,
        "nx": 200,
        "use_control_points_x": False,
        "control_points_x": [
            # Example nonuniform x mesh. To use it set use_uniform_x=False,
            # use_control_points_x=True and keep the first/last point at x0/xmax.
            {"x": 0.0, "dx": 0.010},
            {"x": 0.25, "dx": 0.004},
            {"x": 0.30, "dx": 0.002},
            {"x": 0.40, "dx": 0.004},
            {"x": 1.0, "dx": 0.010},
        ],

        # y direction: 5 cells for the first 2D validation runs
        "use_uniform_y": True,
        "ny": 5,
        "use_control_points_y": False,
        "control_points_y": [
            # Example nonuniform y mesh.
            {"y": 0.0, "dy": 0.0125},
            {"y": 0.025, "dy": 0.00625},
            {"y": 0.05, "dy": 0.0125},
        ],

        "bc_x0": "inlet_left",
        "bc_xmax": "inlet_right",
        "bc_y0": "wall",
        "bc_ymax": "wall",
    },

    "boundary_conditions": {
        "inlet_left": {
            "material_id": 0, "rho": 50.0,
            "u": 0.0, "v": 0.0, "p": 1.0e5,
        },
        "inlet_right": {
            "material_id": 1, "rho": 1000.0,
            "u": 0.0, "v": 0.0, "p": 1.0e9,
        },
        "wall": {"type": "wall"},
    },

    "equations": {
        "materials": {
            0: {"name": "material_0", "eos": "NobleAbelStiffenedGas",
                "nasg": {"gamma": 1.4, "p_c": 0.0, "b": 0.0, "q": 0.0, "cv": 1.0}},
            1: {"name": "material_1", "eos": "NobleAbelStiffenedGas",
                "nasg": {"gamma": 4.4, "p_c": 6.0e8, "b": 0.0, "q": 0.0, "cv": 1.0}},
        },
        "material_transitions": [],
    },

    "initial_condition": {
        "default_state": {
            "material_id": 0, "rho": 50.0,
            "u": 0.0, "v": 0.0, "p": 1.0e5,
        },
        "initial_vapor_region": {
            "enabled": True,
            "x_left": 0.300, "x_right": 1.0,
            "y_bottom": 0.0, "y_top": 0.05,
            "state": {
                "material_id": 1, "rho": 1000.0,
                "u": 0.0, "v": 0.0, "p": 1.0e9,
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
                "frequency": 0,
                "max_iterations": 30,
                "cfl": 0.8,
                "convergence_tolerance": 2.0e-4,
            },
            "reconstruction": {"type": "ThirdOrderUpwindFD", "limiter": "None"},
        }],
    },

    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        # Implemented/tested in 2D now: Extrapolation.
        # "Volume" is kept as a selectable future option, but its existing
        # models/bi_material_state_calculator/volume.py is intentionally untouched.
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
        "type": "Explicit", "max_time": 2.4e-4, "cfl": 0.4,
        "explicit": {"type": "RungeKutta2"},
    },

    "output": {
        "folder_name": "results_rallu_5_2d/",
        "prefix": "results_rallu_5_2d/",
        "solution_name": "solution",
        "time_interval": 2.4e-5,
        "fields": {
            "density": True, "velocity": True, "pressure": True,
            "temperature": False, "sound_speed": False,
            "material_id": True, "levelset0": True,
        },
        "line_plots": [{
            "enabled": True,
            "file_name": "line_center",
            "num_points": 200,
            "x0": 0.0, "xmax": 1.0,
            "y0": 0.025, "ymax": 0.025,
            "time_interval": 2.4e-5,
        }],
        "integrations": [{
            "enabled": True, "every_time_step": True,
            "volume_file": "volume.txt", "mass_file": "mass.txt",
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
