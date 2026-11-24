# 1D expanding high-pressure slab: a four-cell, bubble-like region.
# Non-dimensional units (rho, p, x, t).
# Initially: ambient | 4 high-pressure cells | ambient.
# Both interfaces should move outward with approximately opposite velocities.
# In Cartesian 1D this is a slab, not a spherical bubble.
# This is a qualitative two-interface test (not a conservative reference).

CONFIG = {
    "mesh": {
        "type": "OneDimensional",
        "x0": 0.0,
        "xmax": 1.0,
        "use_uniform": True,
        "nx": 200,                        # dx = 0.01
        "use_control_points": False,
        "control_points_x": [],
        "ny": 1,
        "nz": 1,
        "bc_x0": "inlet",
        "bc_xmax": "inlet2",
    },

    "boundary_conditions": {
        "inlet": {
            "material_id": 0,
            "rho": 1.0, "u": 0.0, "v": 0.0, "w": 0.0, "p": 1.0,
        },
        "inlet2": {
            "material_id": 0,
            "rho": 1.0, "u": 0.0, "v": 0.0, "w": 0.0, "p": 1.0,
        },
    },

    "equations": {
        "materials": {
            0: {
                "name": "ambient_fluid",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4, "p_c": 0.0, "b": 0.0, "q": 0.0,
                    "cv": 1.0,
                },
            },
            1: {
                "name": "dense_pressurized_core",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4, "p_c": 0.0, "b": 0.0, "q": 0.0,
                    "cv": 1.0,
                },
            },
        },
        "material_transitions": [],
    },

    "initial_condition": {
        "default_state": {
            "material_id": 0,
            "rho": 1.0, "u": 0.0, "v": 0.0, "w": 0.0, "p": 1.0,
        },
        # With 100 uniform cells, centers 0.485, 0.495, 0.505, 0.515
        # are the ONLY four cell centers in the core interval below.
        "initial_vapor_region": {
            "enabled": True,
            "x_left": 0.48,
            "x_right": 0.52,
            "state": {
                "material_id": 1,
                "rho": 5.0, "u": 0.0, "v": 0.0, "w": 0.0, "p": 10.0,
            },
        },
    },

    "space": {
        "navier_stokes": {
            "flux": "Rusanov",                # Can also use "Roe"
            "reconstruction": {
                "type": "Constant",             # Can also use "Linear"
                "variable_type": "Conservative",
                "limiter": "GeneralizedMinMod",
                "generalized_minmod_coefficient": 1.2,
            },
        },
        "level_set": [
            {
                "material_id": 1,
                "reinitialization": {"frequency": 0},
                "reconstruction": {"type": "Linear", "limiter": "None"},
            }
        ],
    },

    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        "phase_change": "Volume",  # Two-interface test with nearest-interface velocity.
        "relaxation": {},
    },

    "time": {
        "type": "Explicit",
        "max_time": 0.1,
        "cfl": 0.4,
        "explicit": {"type": "RungeKutta2"},
    },

    "output": {
        "folder_name": "results_expanding_slab/",
        "prefix": "results_expanding_slab/",
        "solution_name": "solution",
        "time_interval": 0.005,
        "fields": {
            "density": True,
            "velocity": True,
            "pressure": True,
            "temperature": False,
            "sound_speed": False,
            "material_id": True,
            "levelset0": True,
        },
        "line_plots": [
            {
                "enabled": True,
                "file_name": "line",
                "num_points": 100,
                "x0": 0.0,
                "xmax": 1.0,
                "y0": 0.0,
                "ymax": 0.0,
                "z0": 0.0,
                "zmax": 0.0,
                "time_interval": 0.025,
            }
        ],
        "integrations": [
            {
                "enabled": True,
                "every_time_step": True,
                "volume_file": "volume.txt",
                "mass_file": "mass.txt",
                "momentum_file": "momentum.txt",
                "total_energy_file": "totalenergy.txt",
                "internal_energy_file": "internalenergy.txt",
            }
        ],
        "log": {
            "enabled": True,
            "file_name": "run.log",
            "log_every_step": False,
            "material_changes": True,
            "material_change_states": False,
            "interface_fluxes": False,
        },
    },
}


def get_config():
    return CONFIG
