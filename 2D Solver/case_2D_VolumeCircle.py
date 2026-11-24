"""Baseline 2D conservative Volume circle test (restored missing example)."""

import copy
from case_2D_CircleAdvection import get_config as circle_config


def get_config():
    cfg = copy.deepcopy(circle_config())
    cfg["mesh"]["nx"] = 32
    cfg["mesh"]["ny"] = 32
    cfg["parallel"]["enabled"] = False
    cfg["multiphase"]["phase_change"] = "Volume"
    cfg["multiphase"]["volume_solver"] = {
        "merge_area_fraction": 0.001,
        "max_hops": 2,
        "max_interface_cfl": 0.5,
        "max_retries": 8,
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
    cfg["time"]["cfl"] = 0.3
    cfg["time"]["max_time"] = 0.02
    cfg["output"]["folder_name"] = "results_volume_circle_2d/"
    cfg["output"]["prefix"] = "results_volume_circle_2d/"
    cfg["output"]["time_interval"] = 0.01
    return cfg
