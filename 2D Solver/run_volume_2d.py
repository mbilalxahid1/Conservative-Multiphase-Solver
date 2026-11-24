"""Run the new integrated Volume case, using the standard 2D driver.

    python run_volume_2d.py
"""

from case_2D_VolumeCircle import get_config
from run_case_2d import run


if __name__ == "__main__":
    result = run(get_config())
    print(f"2D Volume complete: steps={result['step']}, t={result['t']:.8e}")
