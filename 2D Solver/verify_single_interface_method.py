"""Run and record the selectable SingleInterfaceNormal shock-bubble case.

Examples:
  python verify_single_interface_method.py --mesh 72 --target-us 4
  python verify_single_interface_method.py --mesh 320 --target-us 4

The script preserves the supplied case's EOS, CFL, level-set scheme,
interface cut-cell method, small-cell settings and Numba acceleration. Only
mesh, final time, output frequency and file paths are adapted for verification.
It does NOT hide failures or claim reaching target time after failure.
"""
import argparse
import copy
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from case_2D_ShockBubbleCollapse import get_config
from run_case_2d import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, default=72)
    parser.add_argument('--target-us', type=float, default=4.)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--folder', type=Path, default=None)
    args = parser.parse_args()
    if args.mesh < 8 or args.target_us <= 0 or args.threads < 1:
        parser.error('mesh >= 8, target-us > 0, and threads >= 1 required')

    cfg = copy.deepcopy(get_config())
    cfg['mesh']['nx'] = cfg['mesh']['ny'] = args.mesh
    cfg['time']['max_time'] = args.target_us * 1e-6
    cfg['multiphase']['volume_solver']['interface_riemann_method'] = 'SingleInterfaceNormal'
    cfg['parallel'].update({'enabled':True, 'backend':'numba',
                            'threads':args.threads, 'min_cells':0})
    folder = args.folder or Path(f'verification_single_normal_{args.mesh}x{args.mesh}')
    folder = folder.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    cfg['output']['folder_name'] = str(folder)
    cfg['output']['prefix'] = str(folder)
    cfg['output']['time_interval'] = cfg['time']['max_time']
    cfg['output']['line_plots'] = []
    cfg['output']['integrations'] = []
    cfg['output']['log']['enabled'] = False
    cfg['multiphase']['volume_solver']['riemann_fallback_debug'] = False
    info = {'mesh':[args.mesh,args.mesh], 'target_time_s':cfg['time']['max_time'],
            'method':'SingleInterfaceNormal', 'threads':args.threads}
    tic = time.monotonic()
    try:
        result = run(cfg, max_steps=50000, quiet=False)
        info.update(status='completed', time_s=float(result['t']),
                    steps=int(result['step']))
    except Exception as exc:
        info.update(status='failed', exception_type=type(exc).__name__,
                    error=str(exc))
        traceback.print_exc()
    finally:
        info['wall_s'] = time.monotonic()-tic
        path = folder/'volume_conservation.txt'
        if path.exists():
            data = np.loadtxt(path, skiprows=1, ndmin=2)
            if len(data):
                info['accepted_steps'] = int(data[-1,0])
                info['accepted_time_s'] = float(data[-1,1])
                info['max_mass_budget_error'] = float(np.max(np.abs(data[:,11:13])))
                info['max_directional_pressure_difference'] = float(np.max(np.abs(data[:,-1])))
        (folder/'verification_result.json').write_text(json.dumps(info,indent=2))
        print('VERIFICATION_RESULT='+json.dumps(info),flush=True)
    return 0 if info['status']=='completed' else 1


if __name__ == '__main__':
    sys.exit(main())
