#!/usr/bin/env python3
"""
gen_Sources_uniform48_teleseismic.py — Job A/C uniform48 震源（Lin et al. 2014 SKS 準則）
====================================================================================
設計準則（Tier-1 BAZ 掃描，8 × 45° bin 比較）：
  - Epicentral distance Δ = 88° / 100° / 115°（Lin SKS 最佳帶；3 tiers × 16 az）
  - 震源深度 > 50 km（深震，near-field 較弱）：100 / 300 / 600 km
  - 16 azimuths（22.5° 步進）→ 8 bins × 2 directions

輸出至 PSI_D_HFFK/psi_input/（非 project folder）：
  Sources_uniform48.dat
  Sources_uniform48_manifest.csv
  README_uniform48.txt

然後執行：
  python3 scripts/gen_dummy_obs.py

用法：
  cd .../PSI_D_HFFK
  python3 scripts/gen_Sources_uniform48_teleseismic.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

# Taiwan network centroid (same as plot scripts)
TW_LON, TW_LAT = 121.0, 24.0
N_DIR = 16
N_BIN = 8
BIN_WIDTH = 45.0

# Lin et al. 2014a,b — SKS teleseismic window, deep events preferred
DELTA_DEG_TIERS = [88.0, 100.0, 115.0]   # epicentral distance (deg)
FOCAL_DEPTH_KM = [100.0, 300.0, 600.0]  # all > 50 km


def move_great_circle(lon0: float, lat0: float, baz_deg: float, dist_deg: float) -> tuple[float, float]:
    la0, lo0, baz, d = map(math.radians, [lat0, lon0, baz_deg, dist_deg])
    la1 = math.asin(math.sin(la0) * math.cos(d) + math.cos(la0) * math.sin(d) * math.cos(baz))
    lo1 = lo0 + math.atan2(math.sin(baz) * math.sin(d) * math.cos(la0),
                           math.cos(d) - math.sin(la0) * math.sin(la1))
    return math.degrees(lo1), math.degrees(la1)


def haversine_deg(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lo1, la1, lo2, la2 = map(math.radians, [lon1, lat1, lon2, lat2])
    a = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * math.degrees(math.asin(min(1.0, math.sqrt(a))))


def compute_baz(rcv_lon: float, rcv_lat: float, src_lon: float, src_lat: float) -> float:
    rl, sl = map(math.radians, [rcv_lat, src_lat])
    dl = math.radians(src_lon - rcv_lon)
    x = math.sin(dl) * math.cos(sl)
    y = math.cos(rl) * math.sin(sl) - math.sin(rl) * math.cos(sl) * math.cos(dl)
    return math.degrees(math.atan2(x, y)) % 360.0


def bin_index(dir_idx: int) -> int:
    return int((dir_idx * 22.5) // BIN_WIDTH) % N_BIN


def generate_rows() -> list[dict]:
    if len(DELTA_DEG_TIERS) != len(FOCAL_DEPTH_KM):
        raise ValueError("DELTA_DEG_TIERS and FOCAL_DEPTH_KM must have same length (3 tiers)")
    rows: list[dict] = []
    sid = 1
    for tier, (delta, fdep) in enumerate(zip(DELTA_DEG_TIERS, FOCAL_DEPTH_KM)):
        for d in range(N_DIR):
            baz_nom = d * (360.0 / N_DIR)
            lon, lat = move_great_circle(TW_LON, TW_LAT, baz_nom, delta)
            baz_act = compute_baz(TW_LON, TW_LAT, lon, lat)
            b = bin_index(d)
            rows.append(dict(
                src_id=sid, dir_idx=d, bin_idx=b,
                bin_lo=b * 45, bin_hi=(b + 1) * 45,
                baz_nom_deg=round(baz_nom, 2), baz_deg=round(baz_act, 2),
                delta_deg=round(delta, 2),
                delta_check=round(haversine_deg(TW_LON, TW_LAT, lon, lat), 2),
                depth_km=fdep, tier=tier,
                lon=lon, lat=lat,
            ))
            sid += 1
    return rows


def write_sources(path: Path, rows: list[dict]):
    lines = []
    for r in rows:
        dep = -abs(r["depth_km"])
        lines.append(f"{r['src_id']}, {r['lon']:.6f}, {r['lat']:.6f}, {dep:.1f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(path: Path, rows: list[dict]):
    fields = ["src_id", "dir_idx", "bin_idx", "bin_lo", "bin_hi",
              "baz_nom_deg", "baz_deg", "delta_deg", "delta_check",
              "depth_km", "tier", "lon", "lat"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_readme(path: Path, rows: list[dict]):
    text = f"""Sources_uniform48.dat — teleseismic SKS tier-1 catalogue
============================================================
Reference: Lin et al. (2014a,b) — Δ≈88–115°, deep focus (>50 km).

Structure: 48 = 16 azimuths × 3 tiers
  - 16 directions every 22.5° → 8 BAZ bins (45°) × 2 dirs per bin
  - Tier 0: Δ={DELTA_DEG_TIERS[0]}°, depth={FOCAL_DEPTH_KM[0]:.0f} km
  - Tier 1: Δ={DELTA_DEG_TIERS[1]}°, depth={FOCAL_DEPTH_KM[1]:.0f} km
  - Tier 2: Δ={DELTA_DEG_TIERS[2]}°, depth={FOCAL_DEPTH_KM[2]:.0f} km

Centroid: ({TW_LON}°E, {TW_LAT}°N) — Taiwan network reference

Job A/C read:
  ${{PSI_DIR}}/psi_input/Sources_uniform48.dat
  ${{PSI_DIR}}/psi_input/DUMMY_SI_uniform48.dat

After updating sources, rerun:
  sbatch scripts/slurm/jobA2_ray_si_uniform48.sh
  sbatch scripts/slurm/jobC_hffk_uniform48.sh

Sample (src1): Δ={rows[0]['delta_check']}° BAZ={rows[0]['baz_deg']}° bin{rows[0]['bin_idx']}
"""
    path.write_text(text, encoding="utf-8")


def main():
    psi_dir = Path(__file__).resolve().parent.parent
    out_dir = psi_dir / "psi_input"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = generate_rows()
    src_path = out_dir / "Sources_uniform48.dat"
    write_sources(src_path, rows)
    write_manifest(out_dir / "Sources_uniform48_manifest.csv", rows)
    write_readme(out_dir / "README_uniform48.txt", rows)

    # Keep legacy name in sync (deep catalogue = same teleseismic geometry)
    write_sources(out_dir / "Sources_uniform48_deep.dat", rows)

    print(f"PSI_D_HFFK : {psi_dir}")
    print(f"Wrote {src_path}  ({len(rows)} sources)")
    print(f"      {out_dir / 'Sources_uniform48_manifest.csv'}")
    print(f"      {out_dir / 'README_uniform48.txt'}")
    print(f"      {out_dir / 'Sources_uniform48_deep.dat'} (synced)")
    print("\nΔ / depth tiers:")
    for tier, (d, z) in enumerate(zip(DELTA_DEG_TIERS, FOCAL_DEPTH_KM)):
        print(f"  tier{tier}: Δ={d}°  focal depth={z} km  (src {tier*16+1}–{tier*16+16})")
    print("\n8 BAZ bins (2 dirs each):")
    for b in range(N_BIN):
        dirs = [r for r in rows if r["bin_idx"] == b and r["tier"] == 1]
        bazs = [r["baz_deg"] for r in dirs]
        print(f"  bin{b} [{b*45}-{(b+1)*45}°]: BAZ≈{bazs}  Δ=100° tier")
    print("\nNext: python3 scripts/gen_dummy_obs.py")


if __name__ == "__main__":
    main()
