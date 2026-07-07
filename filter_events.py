#!/usr/bin/env python3
"""
filter_events.py — 過濾 PSI_D 觀測資料
條件（Lin et al. 2014a, 2014b）:
  - 震央距離：88° ≤ Δ ≤ 115°（SKS 乾淨視窗）
  - 震源深度：depth > 50 km（減少 near-field 和地殼多重反射）

用法:
  python filter_events.py \
    --sources Sources.dat \
    --receivers Receivers.dat \
    --obs SplittingIntensity_ShearWave.dat \
    --out SplittingIntensity_filtered.dat

輸出: 過濾後的觀測資料（格式與輸入完全相同，可直接餵給 PSI_D）
"""

import argparse
import math
import sys

# ── 大圓距離（Haversine）────────────────────────────────────────────────────
def haversine_deg(lon1, lat1, lon2, lat2):
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
    return math.degrees(2 * math.atan2(math.sqrt(a), math.sqrt(1-a)))

# ── 讀取 Sources.dat ────────────────────────────────────────────────────────
def read_sources(path):
    """id → (lon, lat, elv_km)；elv 為負數表示深度"""
    sources = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            sid  = int(parts[0])
            lon  = float(parts[1])
            lat  = float(parts[2])
            elv  = float(parts[3])          # km；負值 = 深度
            sources[sid] = (lon, lat, elv)
    return sources

# ── 讀取 Receivers.dat ──────────────────────────────────────────────────────
def read_receivers(path):
    """id (string) → (lon, lat, elv_km)"""
    receivers = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            rid = parts[0]
            lon = float(parts[1])
            lat = float(parts[2])
            elv = float(parts[3])
            receivers[rid] = (lon, lat, elv)
    return receivers

# ── 主程式 ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Filter PSI_D observations by distance and depth (Lin 2014)")
    parser.add_argument('--sources',   required=True, help='Sources.dat path')
    parser.add_argument('--receivers', required=True, help='Receivers.dat path')
    parser.add_argument('--obs',       required=True, help='SplittingIntensity .dat path')
    parser.add_argument('--out',       required=True, help='Output filtered .dat path')
    parser.add_argument('--dist_min',  type=float, default=88.0,  help='Min epicentral distance (deg) [default 88]')
    parser.add_argument('--dist_max',  type=float, default=115.0, help='Max epicentral distance (deg) [default 115]')
    parser.add_argument('--depth_min', type=float, default=50.0,  help='Min source depth (km) [default 50]')
    args = parser.parse_args()

    sources   = read_sources(args.sources)
    receivers = read_receivers(args.receivers)

    total = 0
    kept  = 0
    rejected_dist  = 0
    rejected_depth = 0
    rejected_missing = 0

    with open(args.obs) as fin, open(args.out, 'w') as fout:
        for raw_line in fin:
            # 保留空行和注釋（PSI_D 格式不含 header，但以防萬一）
            stripped = raw_line.strip()
            if not stripped or stripped.startswith('#'):
                fout.write(raw_line)
                continue

            total += 1
            parts = [p.strip() for p in stripped.split(',')]

            # PSI_D SplittingIntensity 格式:
            # SI, sigma, period, phase, source_id, receiver_id[, back_azimuth]
            try:
                sid = int(parts[4])
                rid = parts[5]
            except (IndexError, ValueError):
                rejected_missing += 1
                continue

            # 查找 source 和 receiver
            if sid not in sources:
                print(f"  [WARN] source_id {sid} 不在 Sources.dat，跳過", file=sys.stderr)
                rejected_missing += 1
                continue
            if rid not in receivers:
                print(f"  [WARN] receiver_id {rid} 不在 Receivers.dat，跳過", file=sys.stderr)
                rejected_missing += 1
                continue

            src_lon, src_lat, src_elv = sources[sid]
            rec_lon, rec_lat, _       = receivers[rid]

            # 深度（elv 為負數）
            depth_km = abs(src_elv)

            # 震央距離
            dist_deg = haversine_deg(src_lon, src_lat, rec_lon, rec_lat)

            # 過濾條件
            if depth_km <= args.depth_min:
                rejected_depth += 1
                continue
            if dist_deg < args.dist_min or dist_deg > args.dist_max:
                rejected_dist += 1
                continue

            # 通過 → 寫入輸出
            fout.write(raw_line)
            kept += 1

    # ── 統計報告 ──
    print(f"\n{'='*50}")
    print(f"Event-station pair 過濾結果")
    print(f"{'='*50}")
    print(f"  總輸入筆數          : {total}")
    print(f"  距離條件 ({args.dist_min}°–{args.dist_max}°) 刪除 : {rejected_dist}")
    print(f"  深度條件 (>{args.depth_min} km) 刪除    : {rejected_depth}")
    print(f"  找不到 ID 刪除      : {rejected_missing}")
    print(f"  通過（保留）        : {kept}  ({100*kept/total:.1f}%)" if total > 0 else "  (無有效資料)")
    print(f"  輸出檔案            : {args.out}")
    print(f"{'='*50}\n")
    print("引用依據: Lin, Y.-P., Zhao, L., & Hung, S.-H. (2014a). GRL 41, 799–804.")
    print("          Lin, Y.-P., Zhao, L., & Hung, S.-H. (2014b). GRL 41, 8809–8817.")

if __name__ == '__main__':
    main()
