import copy
import numpy as np

from case_2D_ShockBubbleCollapse import get_config
from models.mesh_2d import build_mesh_2d
from models.level_set_solution import initialize_phi_and_material_id_2d
from models.geometry_2d import initial_state_at_point_2d
from models.eos_nasg import primitive_to_conservative_2d, conservative_to_primitive_2d
from models.time_step_cfl import compute_dt_cfl


def main():
    cfg = copy.deepcopy(get_config())
    # Coarse initialization check; physics values and geometry are unchanged.
    cfg["mesh"]["nx"] = 120
    cfg["mesh"]["ny"] = 120
    cfg["parallel"]["enabled"] = False

    x, xf, dx, y, yf, dy = build_mesh_2d(cfg["mesh"])
    phi, material_id = initialize_phi_and_material_id_2d(x, y, dx, dy, cfg)
    U = np.zeros((len(y), len(x), 4), dtype=float)
    mats = cfg["equations"]["materials"]

    for j in range(len(y)):
        for i in range(len(x)):
            s = initial_state_at_point_2d(x[i], y[j], cfg)
            mid = int(s["material_id"])
            U[j, i] = primitive_to_conservative_2d(
                s["rho"], s["u"], s["v"], s["p"], mats[mid])

    # Classify exact initialized primitive states.
    labels = {"postshock_water": 0, "preshock_water": 0, "air_bubble": 0, "other": 0}
    pmin = np.inf; pmax = -np.inf
    for j in range(len(y)):
        for i in range(len(x)):
            mid = int(material_id[j, i])
            rho,u,v,p,*_ = conservative_to_primitive_2d(U[j,i], mats[mid])
            pmin=min(pmin,p); pmax=max(pmax,p)
            if mid == 1:
                labels["air_bubble"] += 1
            elif abs(rho-1323.65) < 1e-8 and abs(u-681.058) < 1e-8 and abs(p-1.9e9) < 1.0:
                labels["postshock_water"] += 1
            elif abs(rho-1000.0) < 1e-8 and abs(u) < 1e-12 and abs(p-1e5) < 1e-3:
                labels["preshock_water"] += 1
            else:
                labels["other"] += 1

    dt = compute_dt_cfl(U, dx, dy, mats, material_id, cfg["time"]["cfl"], cfg["parallel"])
    print("mesh", len(x), len(y))
    print("cell_counts", labels)
    print("material_counts", dict(zip(*np.unique(material_id, return_counts=True))))
    print("phi_min_max", float(phi.min()), float(phi.max()))
    print("pressure_min_max", pmin, pmax)
    print("initial_dt", dt)
    print("shock_sample_left", initial_state_at_point_2d(0.003, 0.012, cfg))
    print("preshock_sample", initial_state_at_point_2d(0.008, 0.012, cfg))
    print("bubble_sample", initial_state_at_point_2d(0.012, 0.012, cfg))

    assert labels["other"] == 0
    assert labels["postshock_water"] > 0
    assert labels["preshock_water"] > 0
    assert labels["air_bubble"] > 0
    assert pmax > 1.0e9
    assert pmin > 0.0
    print("INITIALIZATION_OK")

if __name__ == "__main__":
    main()
