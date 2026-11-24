# 1D two-material shock-tube problem from Arthur Rallu thesis, Fig. 4.15 (Page 135).
# Units: kg, m, s, Pa.
#
# IMPORTANT:
# Only the physical test-case data / final time / output cadence are changed here.
# The solver settings from the supplied Python case are intentionally preserved:
#   - Rusanov flux
#   - Constant reconstruction
#   - Extrapolation interface treatment
#   - RungeKutta2
#   - CFL = 0.4

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

        # Same Python boundary implementation as in the supplied case.
        "bc_x0": "inlet",
        "bc_xmax": "inlet2",
    },

    # -------------------------------------------------------------------------
    # Boundary conditions
    # Rallu Fig. 4.15 / supplied M2C case.
    # -------------------------------------------------------------------------
    "boundary_conditions": {
        "inlet": {
            "material_id": 0,
            "rho": 50.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0e5,
        },
        "inlet2": {
            "material_id": 1,
            "rho": 1000.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0e9,
        },
    },

    # -------------------------------------------------------------------------
    # Equations and materials
    # The Python code uses the NASG form.  Setting b = q = 0 gives the
    # stiffened-gas model used by the supplied M2C case.
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
                    # Required by the current Python EOS dictionary.
                    "cv": 1.0,
                },
            },
            1: {
                "name": "material_1",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 4.4,
                    "p_c": 6.0e8,
                    "b": 0.0,
                    "q": 0.0,
                    # Required by the current Python EOS dictionary.
                    "cv": 1.0,
                },
            },
        },

        "material_transitions": [],
    },

    # -------------------------------------------------------------------------
    # Initial condition
    # M2C Plane[0]: x = 0.3 with normal +x.
    # Therefore:
    #   material 0 occupies x < 0.3,
    #   material 1 occupies x >= 0.3.
    # -------------------------------------------------------------------------
    "initial_condition": {
        "default_state": {
            "material_id": 0,
            "rho": 50.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0e5,
        },

        # The current Python code uses this generic block for material 1.
        # Here material 1 is the high-pressure stiffened-gas/liquid state.
        "initial_vapor_region": {
            "enabled": True,
            "x_left": 0.300,
            "x_right": 1.0,
            "state": {
                "material_id": 1,
                "rho": 1000.0,
                "u": 0.0,
                "v": 0.0,
                "w": 0.0,
                "p": 1.0e9,
            },
        },
    },

    # -------------------------------------------------------------------------
    # Solver settings -- intentionally unchanged from the supplied Python case.
    # -------------------------------------------------------------------------
    "space": {
        "navier_stokes": {
            "flux": "Rusanov",  # Rusanov (LLF) or Roe
            "reconstruction": {
                "type": "Linear",  # Constant or Linear
                "variable_type": "Conservative",  # Primitive, ConservativeCharacteristic, PrimitiveCharacteristic
                "limiter": "GeneralizedMinMod",  # or None
                "generalized_minmod_coefficient": 1.1,
            },
        },

        "level_set": [
            {
                "material_id": 1,
                "reinitialization": {
                    "frequency": 0,
                },
                "reconstruction": {
                    "type": "Linear",
                    "limiter": "None",
                },
            }
        ],
    },

    # -------------------------------------------------------------------------
    # Multiphase solver settings -- intentionally unchanged.
    # -------------------------------------------------------------------------
    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        "phase_change": "Volume",
        "relaxation": {
            # If time_step is omitted, use the number of physical time steps
            # over which the interfacial flux mismatch was accumulated.
            # "time_step": 1000,
        },
    },

    # -------------------------------------------------------------------------
    # Time integration
    # Final time changed to the Fig. 4.15 test value.
    # CFL and RK scheme are intentionally unchanged from your Python case.
    # -------------------------------------------------------------------------
    "time": {
        "type": "Explicit",
        "max_time": 2.0e-4,
        "cfl": 0.4,
        "explicit": {
            "type": "RungeKutta2",
        },
    },

    # -------------------------------------------------------------------------
    # Output
    # Output cadence changed to match the supplied test-case time scale.
    # -------------------------------------------------------------------------
    "output": {
        "folder_name": "results_rallu_5/",
        "prefix": "results_rallu_5/",
        "solution_name": "solution",
        "time_interval": 2.4e-5,

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
                # Write the line solution at the requested final time.
                "time_interval": 2.4e-4,
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
