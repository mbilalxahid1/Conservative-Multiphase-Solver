# Pressure-driven gas bubble expansion in liquid.
# Reference: M2C Tests/ExplosionShyue98/input.st
# Original test: Shyue (1998); M2C input notes also cite Wentao Ma et al., JCP 2023.
#
# Geometry is initialized through models/geometry_2d.py.  This case therefore
# uses the actual circular bubble instead of approximating it by a rectangle.
#
# The test is nondimensional, matching the attached M2C input.st.

CONFIG = {
    # -------------------------------------------------------------------------
    # 2D Cartesian mesh
    # M2C: [0,1] x [0,1], 500 x 500 cells.
    # The original M2C file is formally 3D with one z-cell and symmetry in z;
    # this Python solver is genuinely 2D, so z is omitted.
    # -------------------------------------------------------------------------
    "mesh": {
        "type": "TwoDimensional",
        "x0": 0.0,
        "xmax": 1.0,
        "y0": 0.0,
        "ymax": 1.0,

        "use_uniform_x": True,
        "nx": 200,
        "use_control_points_x": False,
        "control_points_x": [],

        "use_uniform_y": True,
        "ny": 200,
        "use_control_points_y": False,
        "control_points_y": [],

        "bc_x0": "farfield",
        "bc_xmax": "farfield",
        "bc_y0": "farfield",
        "bc_ymax": "farfield",
    },

    # -------------------------------------------------------------------------
    # Far-field state = surrounding liquid state from M2C input.st.
    # material_id is stated explicitly for the Python boundary-condition code.
    # -------------------------------------------------------------------------
    "boundary_conditions": {
        "farfield": {
            "type": "farfield",
            "material_id": 0,
            "rho": 0.991,
            "u": 0.0,
            "v": 0.0,
            "p": 3.059e-4,
        },
    },

    # -------------------------------------------------------------------------
    # Equations of state
    #
    # M2C StiffenedGas:
    #   p = (gamma - 1) rho e - gamma p_c
    # is exactly the b=0, q=0 limit of the Python NASG EOS.
    # cv does not affect this Euler benchmark because temperature is not used.
    # -------------------------------------------------------------------------
    "equations": {
        "materials": {
            0: {
                "name": "liquid",
                "eos": "NobleAbelStiffenedGas",
                "nasg": {
                    "gamma": 5.5,
                    "p_c": 1.505,
                    "b": 0.0,
                    "q": 0.0,
                    "cv": 1.0,
                },
            },
            1: {
                "name": "gas_bubble",
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

    # -------------------------------------------------------------------------
    # Initial condition
    # Background liquid:
    #   rho = 0.991, p = 3.059e-4, u=v=0
    # Circular gas bubble:
    #   center = (0.5,0.5), radius = 0.2
    #   rho = 1.241, p = 2.753, u=v=0
    # -------------------------------------------------------------------------
    "initial_condition": {
        "default_state": {
            "material_id": 0,
            "rho": 0.991,
            "u": 0.0,
            "v": 0.0,
            "p": 3.059e-4,
        },

        "geometry": {
            "enabled": True,
            "type": "Circle",
            "center_x": 0.5,
            "center_y": 0.5,
            "radius": 0.2,
            "state": {
                "material_id": 1,
                "rho": 1.241,
                "u": 0.0,
                "v": 0.0,
                "p": 2.753,
            },
        },
    },

    # -------------------------------------------------------------------------
    # Spatial discretization
    # M2C input.st:
    #   Flux = LocalLaxFriedrichs
    #   Reconstruction = Constant
    # In the current Python code the matching same-material flux is Rusanov.
    # -------------------------------------------------------------------------
    "space": {
        "navier_stokes": {
            "flux": "Rusanov",
            "reconstruction": {
                "type": "Constant",
            },
        },

        "level_set": [{
            "material_id": 1,

            # M2C input.st does not request periodic reinitialization for this case.
            "reinitialization": {
                "frequency": 0,
                "max_iterations": 30,
                "cfl": 0.4,
                "convergence_tolerance": 2.0e-4,
            },

            # Current Python level-set implementation is the third-order
            # directional upwind finite-difference method.
            "reconstruction": {
                "type": "ThirdOrderUpwindFD",
                "limiter": "None",
            },
        }],
    },

    # -------------------------------------------------------------------------
    # Multiphase interface treatment
    # Direct translation of M2C:
    #   Flux = Numerical
    #   ReconstructionAtInterface = Constant
    #   PhaseChange = Extrapolation
    # -------------------------------------------------------------------------
    "multiphase": {
        "flux": "Numerical",
        "reconstruction_at_interface": "Constant",
        "phase_change": "Volume",
        "phase_change_direction": "Upwind",
    },

    # -------------------------------------------------------------------------
    # Time integration
    # M2C uses RK3, MaxTime=0.058, CFL=0.4.
    # Current Python driver accepts SSP-RK2 only, so RK2 is used here while
    # preserving the M2C final time and CFL.
    # -------------------------------------------------------------------------
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
        "max_time": 0.058,
        "cfl": 0.4,
        "explicit": {
            "type": "RungeKutta2",
        },
    },

    # -------------------------------------------------------------------------
    # Output
    # Directly follows the M2C output fields and centerline request.
    # -------------------------------------------------------------------------
    "output": {
        "folder_name": "results_explosion_shyue98_2d/",
        "prefix": "results_explosion_shyue98_2d/",
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

        # M2C LinePlot[0]: x from 0 to 1 at y=0.5, 500 points.
        "line_plots": [{
            "enabled": True,
            "file_name": "line",
            "num_points": 500,
            "x0": 0.0,
            "xmax": 1.0,
            "y0": 0.5,
            "ymax": 0.5,
            "time_interval": 0.01,
        }],

        # Extra Python diagnostics useful for benchmarking mass conservation.
        # These are not part of the original M2C input.st physics.
        "integrations": [{
            "enabled": True,
            "every_time_step": False,
            "time_interval": 0.001,
            "volume_file": "volume.txt",
            "mass_file": "mass.txt",
            "momentum_file": "momentum.txt",
            "total_energy_file": "totalenergy.txt",
            "internal_energy_file": "internalenergy.txt",
        }],

        # Avoid very large per-step log files for the 500 x 500 benchmark.
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
