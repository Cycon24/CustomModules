import argparse, yaml, subprocess, sys
from pathlib import Path

def run(cmd, cwd=None):
    p = subprocess.Popen(cmd, cwd=cwd, shell=False)
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

def main():
    ap = argparse.ArgumentParser(description="Batch 3-D post-processing driver for a single sweep-parameter folder")
    ap.add_argument("--config", required=True, help="Path to post3d_config.yaml")
    ap.add_argument("--param_dir", default=".", help="Path to a single parameter directory containing entire_surface_restart.csv")
    ap.add_argument("--units", required=False, help="Optional plot_units.yaml to override")
    args = ap.parse_args()

    param_dir = Path(args.param_dir).resolve()
    cfg_path = Path(args.config).resolve()
    units = args.units

    run([sys.executable, str(Path(__file__).parent / "extract_probe_surfaces.py"), "--config", str(cfg_path), "--param_dir", str(param_dir)])

    cmd = [sys.executable, str(Path(__file__).parent / "plot_cross_section_contours.py"), "--config", str(cfg_path), "--param_dir", str(param_dir)]
    if units: cmd += ["--units", units]
    run(cmd)

    cmd = [sys.executable, str(Path(__file__).parent / "plot_radial_averages.py"), "--config", str(cfg_path), "--param_dir", str(param_dir)]
    if units: cmd += ["--units", units]
    run(cmd)

if __name__ == "__main__":
    main()
