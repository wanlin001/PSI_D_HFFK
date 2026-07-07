#!/usr/bin/env python3
"""
Batch runner: iterate over all projects/<name>/ dirs and run PSI_D.
Usage:
  python scripts/batch_run.py                  # run all projects
  python scripts/batch_run.py 0608 0609        # run only matching names
"""

import subprocess
import sys
import glob
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JULIA = "julia"

def run_project(proj_dir: Path) -> bool:
    print(f"\n{'='*60}")
    print(f"Project: {proj_dir.name}")
    print(f"{'='*60}")
    result = subprocess.run(
        [JULIA, f"--project={REPO_ROOT}", str(REPO_ROOT / "scripts/run_psi.jl"), str(proj_dir)],
        cwd=REPO_ROOT,
    )
    if result.returncode != 0:
        print(f"  FAILED: {proj_dir.name}")
        return False
    return True

def main():
    filters = sys.argv[1:]
    all_projects = sorted(Path(REPO_ROOT / "projects").glob("*/"))
    projects = [p for p in all_projects if p.is_dir() and (p / "psi_config.toml").exists()]

    if filters:
        projects = [p for p in projects if any(f in p.name for f in filters)]

    if not projects:
        print("No matching projects found.")
        return

    print(f"Running {len(projects)} project(s):")
    for p in projects:
        print(f"  {p.name}")

    failed = []
    for proj in projects:
        ok = run_project(proj)
        if not ok:
            failed.append(proj.name)

    print(f"\n{'='*60}")
    print(f"Done: {len(projects) - len(failed)}/{len(projects)} succeeded")
    if failed:
        print(f"Failed: {', '.join(failed)}")

if __name__ == "__main__":
    main()
