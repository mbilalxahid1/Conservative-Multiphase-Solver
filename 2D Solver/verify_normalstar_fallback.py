"""Standalone short coupled verification for the opt-in normal-star flux fallback."""
import argparse
from pathlib import Path
import numpy as np
from case_2D_ShockBubbleCollapse import get_config
from run_case_2d import run

p=argparse.ArgumentParser()
p.add_argument('--mode',choices=['Riemann','Interface'],default='Riemann')
p.add_argument('--size',type=int,default=72)
p.add_argument('--end-time',type=float,default=2.0e-6)
a=p.parse_args()

cfg=get_config()
cfg['mesh']['nx']=cfg['mesh']['ny']=a.size
cfg['time']['max_time']=a.end_time
cfg['output']['folder_name']=str(Path('verification_shockbubble')/(a.mode+'_'+str(a.size)))
cfg['output']['time_interval']=a.end_time
cfg['parallel']['threads']=2
cfg['multiphase']['volume_solver']['normal_mode']=a.mode
cfg['multiphase']['volume_solver']['normal_star_flux_fallback']=True
cfg['multiphase']['volume_solver']['riemann_fallback_debug']=False
r=run(cfg,max_steps=10000,quiet=True)
errors=np.array([np.max(abs(k['mass_error'])) for k in r['volume_conservation']])
all_budget=np.array([np.max(abs(k['conservation_error'])) for k in r['volume_conservation'] if 'conservation_error' in k])
print(f"RESULT mode={a.mode} cells={a.size} step={r['step']} t={r['t']:.17e} max_mass_error={np.max(errors):.17e} max_budget_error={np.max(all_budget) if len(all_budget) else np.nan:.17e}",flush=True)
