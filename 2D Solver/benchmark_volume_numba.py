"""Compare the original serial and Numba 2D Volume steps on the same state.

Usage: python benchmark_volume_numba.py --size 120 --threads 4
First-use compilation occurs on a small warm-up problem and is reported
separately. No simulation result files are changed or overwritten.
"""

import argparse
import copy
import time

import numpy as np

from case_2D_CircleAdvection import get_config
from models.bi_material_state_calculator.volume import (
    initialize_volume_state_2d, reconstruct_volume_geometry_2d,
    volume_step_2d,
)
from models.eos_nasg import primitive_to_conservative_2d
from models.parallel_numba import configure_numba, NUMBA_AVAILABLE


def make_case(size):
    cfg = copy.deepcopy(get_config())
    cfg["multiphase"]["phase_change"] = "Volume"
    cfg["space"]["level_set"][0]["reinitialization"]["frequency"] = 0
    cfg["parallel"]["min_cells"] = 0
    f = np.linspace(0.0, 1.0, size + 1)
    x = 0.5 * (f[1:] + f[:-1])
    X, Y = np.meshgrid(x, x)
    phi = np.hypot(X - 0.5, Y - 0.5) - 0.23
    ids = (phi < 0).astype(int)
    materials = cfg["equations"]["materials"]
    U = np.empty((size, size, 4))
    for m, density in enumerate((1.0, 2.0)):
        U[ids == m] = primitive_to_conservative_2d(
            density, 0.1, 0.05, 1.0, materials[m])
    geom = reconstruct_volume_geometry_2d(phi, f, f, cfg)
    q0, q1 = initialize_volume_state_2d(U, ids, geom, materials, x, x)
    return cfg, f, phi, q0, q1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, default=120)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if not NUMBA_AVAILABLE:
        raise RuntimeError("Numba is unavailable. Install it before benchmarking.")
    if args.size < 16 or args.threads < 1:
        raise ValueError("Use --size >= 16 and --threads >= 1.")

    cfg, f, phi, q0, q1 = make_case(32)
    cfg["parallel"].update({"enabled": True, "threads": args.threads})
    configure_numba(cfg)
    t0 = time.perf_counter()
    volume_step_2d(q0, q1, phi, f, f, cfg, 5.e-5)
    print(f"Numba warm-up (includes JIT compilation): {time.perf_counter()-t0:.3f} s")

    cfg, f, phi, q0, q1 = make_case(args.size)
    dt = 6.e-5
    timings = {}
    results = {}
    for enabled in (False, True):
        cfg["parallel"]["enabled"] = enabled
        t0 = time.perf_counter()
        result = volume_step_2d(q0, q1, phi, f, f, cfg, dt)
        elapsed = time.perf_counter() - t0
        label = "Numba" if enabled else "Serial"
        results[label] = result
        timings[label] = elapsed
        print(f"{label:>6}: {elapsed:8.3f} s, "
              f"material-mass error={np.max(abs(result[4]['mass_error'])):.3e}")
    print(f"Speedup: {timings['Serial']/timings['Numba']:.1f}x")
    for index, label in enumerate(("Q0", "Q1", "phi")):
        a, b = results["Serial"][index], results["Numba"][index]
        print(f"Maximum absolute {label} difference: {np.max(abs(a-b)):.3e}")


if __name__ == "__main__":
    main()
