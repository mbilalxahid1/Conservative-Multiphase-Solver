# 1D two-material shock-tube problem from Rallu thesis, Fig. 4.12.
# Non-dimensional units.

CONFIG = {
    # -------------------------------------------------------------------------
    # Mesh
    # -------------------------------------------------------------------------
    "mesh": {
        "type": "OneDimensional",
        "x0": 0.0,
        "xmax": 1.0,

        "use_uniform": True,
        "nx": 200,

        "use_control_points": False,
        "control_points_x": [],

        "ny": 1,
        "nz": 1,

        # Temporary choice for the current Python boundary implementation.
        "bc_x0": "inlet",
        "bc_xmax": "inlet2",
    },

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    "boundary_conditions": {
        "inlet": {
            "material_id": 0,
            "rho": 1.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0,
        },
        "inlet2": {
            "material_id": 1, # material ID 1 is called as vapor in code
            "rho": 0.125,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 0.1,
        },
    },

    # -------------------------------------------------------------------------
    # Equations and materials
    # p_c = b = q = 0 reduces the NASG model to the perfect-gas form used here.
    # cv is required by the current EOS dictionary but does not affect rho, u, p.
    # -------------------------------------------------------------------------
    "equations": {
        "materials": {
            0: {
                "name": "material_0",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.4,
                    "p_c": 0.0,
                    "b": 0.0,
                    "q": 0.0,
                    "cv": 1.0,
                },
            },
            1: {
                "name": "material_1",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.2,
                    "p_c": 0.0,
                    "b": 0.0,
                    "q": 0.0,
                    "cv": 1.0,
                },
            },
        },

        # No thermally driven phase transition in this shock-tube problem.
        "material_transitions": [],
    },

    # -------------------------------------------------------------------------
    # Initial condition
    # Material 0 occupies x < 0.5.
    # Material 1 occupies x >= 0.5.
    # -------------------------------------------------------------------------
    "initial_condition": {
        "default_state": {
            "material_id": 0,
            "rho": 1.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0,
        },

        # The current code uses this block as the generic material-1 region.
        "initial_vapor_region": {
            "enabled": True,
            "x_left": 0.500,
            "x_right": 1.0,
            "state": {
                "material_id": 1,
                "rho": 0.125,
                "u": 0.0,
                "v": 0.0,
                "w": 0.0,
                "p": 0.1,
            },
        },
    },

    # -------------------------------------------------------------------------
    # Space discretization
    # The current driver actually uses constant reconstruction and Rusanov flux.
    # -------------------------------------------------------------------------
    "space": {
        "navier_stokes": {
            "flux": "Rusanov",
            "reconstruction": {
                "type": "Constant",
            },
        },

        "level_set": [
            {
                "material_id": 1,
                "reinitialization": {
                    # Disable periodic reinitialization for this test.
                    "frequency": 0,
                },
                # This block is retained for compatibility with the case format.
                "reconstruction": {
                    "type": "Linear",
                    "limiter": "None",
                },
            }
        ],
    },

    # -------------------------------------------------------------------------
    # Multiphase interface treatment
    # The current code uses its bimaterial interface flux and donor extrapolation.
    # -------------------------------------------------------------------------
    "multiphase": {
    "flux": "Numerical",
    "reconstruction_at_interface": "Constant",
    "phase_change": "Volume",
    "relaxation": {
        # If time_step is omitted, use the number of physical time steps
        # over which the interfacial flux mismatch was accumulated.
        #"time_step": 1000,
        },
    },

    # -------------------------------------------------------------------------
    # Time integration
    # -------------------------------------------------------------------------
    "time": {
        "type": "Explicit",
        "max_time": 0.15,
        "cfl": 0.4,
        "explicit": {
            "type": "RungeKutta2",
        },
    },

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------
    "output": {
        "folder_name": "results_rallu/",
        "prefix": "results_rallu/",
        "solution_name": "solution",
        "time_interval": 0.01,

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
                "num_points": 200,
                "x0": 0.0,
                "xmax": 1.0,
                "y0": 0.0,
                "ymax": 0.0,
                "z0": 0.0,
                "zmax": 0.0,
                "time_interval": 0.15,
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
            "log_every_step": True,
            "material_changes": True,
            "material_change_states": True,
            "interface_fluxes": False,
        },
    },
}


def get_config():
    return CONFIG
