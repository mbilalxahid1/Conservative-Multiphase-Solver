# Case configuration
CONFIG = {
    # Mesh
    "mesh": {
        "type": "OneDimensional",
        "x0": -20.0,     # [mm]
        "xmax": 20.0,    # [mm]

        # --- Uniform mesh ---
        "use_uniform": True, # or False for non-uniform
        "nx": 1000,      # NumberOfCellsX (uniform)


        # --- Nonuniform mesh via control points for later use if I need
        "use_control_points": False,
        "control_points_x": [
            # for each control point: {"x": Coordinate [mm], "dx": CellWidth [mm]}
            {"x": -20.0, "dx": 0.1},
            {"x":  -8.0, "dx": 0.05},
            {"x":  -4.0, "dx": 0.025},
            {"x":  -2.5, "dx": 0.015},
            {"x":   0.0, "dx": 0.015},
            {"x":   2.5, "dx": 0.015},
            {"x":   4.0, "dx": 0.025},
            {"x":   8.0, "dx": 0.05},
            {"x":  20.0, "dx": 0.1},
        ],
        # Always keep these as 1 
        "ny": 1,
        "nz": 1,
        
        # Boundary types
        "bc_x0": "inlet",
        "bc_xmax": "inlet",
    },
    
# -------------------------------------------------------------------------
# Boundary conditions (for outlet and wall type B.Cs, we do not need to specify any values)
# -------------------------------------------------------------------------
    "boundary_conditions": {
        "inlet": {
            "rho": 0.0010,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "p": 1.0e5,
        }
    },
    
    
    # Governing equations / Materials / Phase transition
    "equations": {
        "materials": {
            0: {  # water (liquid)
                "name": "water_liquid",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 2.057,           # values from Zhao JFM 2024 Paper
                    "p_c": 1.066e9,           # [Pa]
                    "b": 0.0,                 # [1/rho units consistent with solver]
                    "q": -1994.674e9,         # [mm^2/s^2]
                    "cv": 3.449e9,            # [mm^2/(s^2*K)]
                    # "rho_cutoff": 1.0e-12,  # if need to stop crashing
                    # "p_cutoff": -3.0e8,
                },
            },

            1: {  # water (vapor)
                "name": "water_vapor",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 1.327,
                    "p_c": 0.0,               # [Pa]
                    "b": 0.0,
                    "q": 1995.0e9,            # [mm^2/s^2]
                    "cv": 1.2e9,              # [mm^2/(s^2*K)]
                    # "rho_cutoff": 1.0e-13,
                    # "p_cutoff": 1.0e-13,
                },
            },
        },

        "material_transitions": [   # only vaporization, did not used for simple shock cases
            {   # vaporization
                "name": "vaporization",
                "from_material_id": 0,
                "to_material_id": 1,
                "T_upper": 373.15,        # [K] assumed vaporization temperature
                "latent_heat": 2.2564e12, # [mm^2/s^2] (as given)

                "kinetics_model": "None",
                "constant_kinetics": {
                    "duration": 1.0e-13,  # [s] instantaneous
                },
            }
        ],
    },

    # -------------------------------------------------------------------------
    # Initial condition
    # -------------------------------------------------------------------------
    "initial_condition": {
    "default_state": {
        "material_id": 0,
        "rho": 0.001,      # [g/mm^3]
        "u": 0.0,          # [mm/s] (or consistent with solver units)
        "v": 0.0,
        "w": 0.0,
        "p": 1.0e5,        # [Pa]
        # Optional: may later add "T": ...
    },

    "initial_vapor_region": {
        "enabled": True, 
        # 1D interval [x_left, x_right]
        "x_left": -0.5,     # [mm]
        "x_right": 1.5,     # [mm]

        # Vapor state inside initial vapor
        "state": {
            "material_id": 1,
            "rho": 1.0e-5,   # [g/mm^3]
            "u": 0.0,        # [mm/s]
            "v": 0.0,
            "w": 0.0,
            "p": 1.0e8,      # [Pa], High pressure region for bubble expansion
        },
    },

    # Optional geometry-based initialization (COMMENTED) # Not worked on this # For future
    # "geometric_entities": {
    #     "cylinder_and_cone": [
    #         {
    #             "axis": (1.0, 0.0, 0.0),
    #             "base_center": (-0.2, 0.0, 0.0),
    #             "cylinder_radius": 0.1,
    #             "cylinder_height": 0.4,
    #             "cone_height": 0.0,
    #             "side": "Interior",
    #             "initial_state": {
    #                 "material_id": 0,
    #                 "rho": 0.001,
    #                 "u": 0.0, "v": 0.0, "w": 0.0,
    #                 "p": 1.0e5,
    #             },
    #         }
    #     ]
    # },
},


    # -------------------------------------------------------------------------
    # Space discretization / numerics
    # -------------------------------------------------------------------------
    "space": {
        "navier_stokes": {
            "flux": "LLF", # only type available
            "reconstruction": {
                "type": "Linear",
                "variable_type": "ConservativeCharacteristic",  # C-C = C = P-C > P
                "slope_near_interface": "NonZero",
                "limiter": "GeneralizedMinMod",
                "generalized_minmod_coefficient": 1.0,
            },
        },

        # Level-set settings for the current 1D interface tracker
        # What is actually used right now: #other are just dummy for now
        # - material_id
        # - reinitialization["frequency"]
        "level_set": [
            {
                "material_id": 1,
                "bandwidth": 6,
                "reinitialization": {
                    "frequency": 2,
                    "max_its": 800,
                    "first_layer_treatment": "Fixed",
                    "cfl": 0.4,
                    "convergence_tolerance": 1.0e-5,
                },
                "reconstruction": {
                    "type": "Linear",
                    "limiter": "None",
                },
            }
        ],
    },

    # -------------------------------------------------------------------------
    # MultiPhase (COMMENTED for now)
    # -------------------------------------------------------------------------
    "multiphase": {
    "flux": "Numerical",
    "reconstruction_at_interface": "Constant",
    "phase_change": "Volume",
    "relaxation": {
        # If time_step is omitted, use the number of physical time steps
        # over which the interfacial flux mismatch was accumulated.
        #"time_step": 10,
        },
    },

    # -------------------------------------------------------------------------
    # Time integration
    # -------------------------------------------------------------------------
    "time": {
        "type": "Explicit",
        "max_time": 5.0e-6,  # [s]
        "cfl": 0.4,
        "explicit": {
            "type": "RungeKutta2",
        },
        # "max_its": 500,
    },

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------
    "output": {
        "folder_name": "results/",
        "solution_name": "solution",
        "time_interval": 1.0e-7,

        "fields": {
    "density": True,
    "velocity": True,
    "pressure": True,
    "temperature": True,
    "sound_speed": True,
    "material_id": True,
    "levelset0": True,
    },

        "line_plots": [
            {
                "enabled": True,
                "file_name": "line",
                "num_points": 1000,
                "x0": -20.0,
                "xmax": 20.0,
                "y0": 0.0,
                "ymax": 0.0,
                "z0": 0.0,
                "zmax": 0.0,
                #"time_interval": 0.15,
            }
        ],

        "integrations": [
            {
                "enabled": True,
                "every_time_step": False, # Set every_time_step=True to write at every solver time step.
                "time_interval": 1.0e-7,
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
            "material_changes": True,
            "material_change_states": True,
            "interface_fluxes": True,
        },
    },
}


def get_config():
    return CONFIG