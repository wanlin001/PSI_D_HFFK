"""
SinkingSlab HFFK Period Effect — 3D Spatial Visualization
==========================================================
用法：
  python plot_si_3d.py \
    --ray   psi_output/SinkingSlab_ray/SYN_SplittingIntensity_ShearWave.dat \
    --hffk3 psi_output/SinkingSlab_hffk3s/SYN_SplittingIntensity_ShearWave.dat \
    --hffk8 psi_output/SinkingSlab_hffk8s/SYN_SplittingIntensity_ShearWave.dat \
    --hffk25 psi_output/SinkingSlab_hffk25s/SYN_SplittingIntensity_ShearWave.dat \
    [--receivers psi_input/Receivers.dat]  # 站點座標（default: auto-detect）
    [--src 5]           # 只畫單一震源（default: 全部）
    [--mode map|3d|surface|diff|bar|all]
    [--vmax 0.4]        # colorbar 上下限
    [--out figs/prefix] # 存圖前綴（省略 → interactive）

PSI_D SplittingIntensity 輸出欄位：
  col0: SI_syn  col1: SI_obs  col2: ?  col3: phase  col4: evid
  col5: station_name  col6: ???  col7: azimuth(rad)
站點座標另存於 Receivers.dat：name, lon, lat, elev
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401
from pathlib import Path


# ── helpers ──────────────────────────────────────────────────────────────────

def load_receivers(path):
    """讀 Receivers.dat → {station_name: (lon, lat)}"""
    rec = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = [x.strip() for x in ln.split(',')]
            rec[parts[0]] = (float(parts[1]), float(parts[2]))
    return rec


def find_receivers(si_path):
    """從 SI 輸出路徑往上找 psi_input/Receivers.dat"""
    p = Path(si_path).resolve()
    for parent in p.parents:
        candidate = parent / 'psi_input' / 'Receivers.dat'
        if candidate.exists():
            return str(candidate)
    return None


def load_sp(path, src=None, receivers=None):
    """
    讀 PSI_D SplittingParameters 輸出，join Receivers.dat 取座標。
    輸出欄位：dt_syn, phi_syn(rad), dt_err, phi_err, weight, phase, evid, station, channel, az
    回傳 dict: stlon, stlat, dt, phi_rad
    """
    rows = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = [x.strip() for x in ln.split(',')]
            dt      = float(parts[0])
            phi_rad = float(parts[1])
            evid    = int(parts[6])
            station = parts[7]
            rows.append((dt, phi_rad, evid, station))

    if src is not None:
        rows = [r for r in rows if r[2] == src]

    # group by station: circular mean of phi, mean of dt
    from collections import defaultdict
    sta_data = defaultdict(lambda: {'dt': [], 'phi': []})
    for dt, phi_rad, evid, sta in rows:
        sta_data[sta]['dt'].append(dt)
        sta_data[sta]['phi'].append(phi_rad)

    stations = list(sta_data.keys())
    dt_arr  = np.array([np.mean(sta_data[s]['dt']) for s in stations])
    # circular mean: average unit vectors, then atan2
    phi_mean = np.array([
        np.arctan2(np.mean(np.sin(sta_data[s]['phi'])),
                   np.mean(np.cos(sta_data[s]['phi'])))
        for s in stations
    ])
    # circular std: sqrt(-2 * ln(R)) where R = mean resultant length
    phi_std = np.array([
        np.sqrt(-2 * np.log(np.clip(
            np.hypot(np.mean(np.sin(sta_data[s]['phi'])),
                     np.mean(np.cos(sta_data[s]['phi']))), 1e-9, 1.0)))
        for s in stations
    ])

    stlon = np.full(len(stations), np.nan)
    stlat = np.full(len(stations), np.nan)
    if receivers:
        for i, sta in enumerate(stations):
            if sta in receivers:
                stlon[i], stlat[i] = receivers[sta]

    return {'stlon': stlon, 'stlat': stlat, 'dt': dt_arr, 'phi': phi_mean, 'phi_std': phi_std}


def load_si(path, src=None, receivers=None):
    """
    讀 PSI_D SplittingIntensity 輸出，join Receivers.dat 取座標，回傳 dict of arrays。
    輸出欄位：SI_syn, SI_obs, ?, phase, evid, station, ???, azimuth
    """
    rows = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = [x.strip() for x in ln.split(',')]
            si_val  = float(parts[0])
            evid    = int(parts[4])
            station = parts[5]
            rows.append((si_val, evid, station))

    if src is not None:
        rows = [(s, e, st) for s, e, st in rows if e == src]

    si_arr   = np.array([r[0] for r in rows])
    evid_arr = np.array([r[1] for r in rows])
    sta_arr  = [r[2] for r in rows]

    # join coordinates
    stlon = np.full(len(rows), np.nan)
    stlat = np.full(len(rows), np.nan)
    if receivers:
        for i, sta in enumerate(sta_arr):
            if sta in receivers:
                stlon[i], stlat[i] = receivers[sta]

    return {
        'evid':    evid_arr,
        'station': sta_arr,
        'stlon':   stlon,
        'stlat':   stlat,
        'si':      si_arr,
    }


def diverging_cmap():
    return plt.cm.RdBu_r


def make_norm(si, vmax=None):
    if vmax is None:
        vmax = np.percentile(np.abs(si), 98)
        vmax = max(vmax, 1e-6)
    return mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)


# ── plot functions ────────────────────────────────────────────────────────────

def plot_2d_maps(datasets, titles, sp_datasets=None, fig_kw=None,
                 vmax=None, shared_scale=False, quiver_skip=3,
                 quiver_fixed_len=False, out=None):
    """
    2D 散點圖（顏色=SI）+ 可選 quiver（快軸方向箭頭，來自 SplittingParameters）。
    sp_datasets: list of dict from load_sp()，長度須與 datasets 相同（None 代表該子圖不畫箭頭）
    quiver_skip: 每隔幾個站點畫一支箭頭（避免太密）
    """
    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(4.5*n, 4.8), **(fig_kw or {}))
    if n == 1:
        axes = [axes]
    cmap = diverging_cmap()

    if shared_scale:
        all_si = np.concatenate([d['si'] for d in datasets])
        norms = [make_norm(all_si, vmax)] * n
    else:
        norms = [make_norm(d['si'], vmax) for d in datasets]

    for i, (ax, d, title, norm) in enumerate(zip(axes, datasets, titles, norms)):
        mask = ~np.isnan(d['stlon'])
        sc = ax.scatter(d['stlon'][mask], d['stlat'][mask], c=d['si'][mask],
                        cmap=cmap, norm=norm, s=22, edgecolors='k', linewidths=0.2)

        # quiver：快軸方向箭頭（若有 SP 資料）
        if sp_datasets is not None and i < len(sp_datasets) and sp_datasets[i] is not None:
            sp = sp_datasets[i]
            sp_mask = ~np.isnan(sp['stlon'])
            idx = np.where(sp_mask)[0][::quiver_skip]
            phi = sp['phi'][idx]
            dt  = sp['dt'][idx]
            # 畫線段（tick mark 風格），完全控制長度
            if quiver_fixed_len:
                half_len = np.full_like(dt, 0.15)
            else:
                all_dt = np.concatenate([s['dt'] for s in sp_datasets if s is not None])
                dt_max = np.percentile(all_dt, 95)
                half_len = np.clip(dt / dt_max, 0.20, 1.0) * 0.15
            dx = half_len * np.cos(phi)
            dy = half_len * np.sin(phi)
            phi_std = sp['phi_std'][idx] if 'phi_std' in sp else np.zeros_like(phi)
            for xi, yi, r, dxi, dyi, phi_i, std_i in zip(
                    sp['stlon'][idx], sp['stlat'][idx], half_len,
                    dx, dy, phi, phi_std):
                # 1-std fan on both ends (180° ambiguity symmetric)
                for sign in (1, -1):
                    fan_a = np.linspace(phi_i + sign*std_i, phi_i - sign*std_i, 20)
                    fx = np.concatenate([[xi], xi + sign * r * np.cos(fan_a)])
                    fy = np.concatenate([[yi], yi + sign * r * np.sin(fan_a)])
                    ax.fill(fx, fy, color='k', alpha=0.12, lw=0)
                ax.plot([xi - dxi, xi + dxi], [yi - dyi, yi + dyi],
                        'k-', lw=0.5, solid_capstyle='butt')

        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Lon (°)')
        ax.set_ylabel('Lat (°)')
        ax.set_aspect('equal', 'box')
        std_val = np.std(d['si'])
        ax.text(0.02, 0.98,
                f'std={std_val:.4f}\n|SI|≤{norm.vmax:.4f}',
                transform=ax.transAxes, va='top', fontsize=7.5,
                bbox=dict(fc='white', alpha=0.75, ec='none'))
        fig.colorbar(sc, ax=ax, label='SI', shrink=0.75, pad=0.02)

    title_suffix = ' + fast axis (quiver)' if sp_datasets else ''
    fig.suptitle(f'SinkingSlab: Ray Theory vs HFFK{title_suffix}\n(independent color scale per panel)',
                 fontsize=11, y=1.02)
    plt.tight_layout()
    _save_or_show(fig, out)


def plot_3d_scatter(datasets, titles, vmax=None, out=None):
    """
    3D 散點圖：x=stlon, y=stlat, z=SI，顏色=SI。
    每個 dataset 一個子圖。
    """
    n = len(datasets)
    fig = plt.figure(figsize=(5*n, 5))
    cmap = diverging_cmap()
    norms = [make_norm(d['si'], vmax) for d in datasets]  # per-panel scale

    for i, (d, title, norm) in enumerate(zip(datasets, titles, norms)):
        ax = fig.add_subplot(1, n, i+1, projection='3d')
        colors = cmap(norm(d['si']))
        ax.scatter(d['stlon'], d['stlat'], d['si'],
                   c=colors, s=20, depthshade=False)
        ax.set_xlabel('Lon', fontsize=8, labelpad=2)
        ax.set_ylabel('Lat', fontsize=8, labelpad=2)
        ax.set_zlabel('SI', fontsize=8, labelpad=2)
        ax.set_title(title, fontsize=10)
        ax.tick_params(labelsize=7)
        # 畫 z=0 平面（參考線）
        lon_ok = d['stlon'][~np.isnan(d['stlon'])]
        lat_ok = d['stlat'][~np.isnan(d['stlat'])]
        if len(lon_ok) and len(lat_ok):
            llon = [lon_ok.min(), lon_ok.max()]
            llat = [lat_ok.min(), lat_ok.max()]
            xx, yy = np.meshgrid(llon, llat)
            ax.plot_surface(xx, yy, np.zeros_like(xx), alpha=0.08, color='gray')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='SI', shrink=0.6, pad=0.05)
        ax.text2D(0.02, 0.98, f'std={d["si"].std():.4f}',
                  transform=ax.transAxes, fontsize=7.5, va='top',
                  bbox=dict(fc='white', alpha=0.7, ec='none'))

    fig.suptitle('SinkingSlab 3D: station position vs SI (per-panel scale)', fontsize=12)
    plt.tight_layout()
    _save_or_show(fig, out)


def plot_3d_surface(d_ray, d_hffk, title_hffk='HFFK', vmax=None, out=None):
    """
    3D surface/trisurf：把 (stlon, stlat, SI) 做成曲面，
    左：ray theory，右：HFFK。
    需要站點在規則格點上（SinkingSlab 合成資料符合此條件）。
    """
    from matplotlib.tri import Triangulation

    fig = plt.figure(figsize=(12, 5))
    cmap = diverging_cmap()
    all_si = np.concatenate([d_ray['si'], d_hffk['si']])
    norm = make_norm(all_si, vmax)

    for i, (d, title) in enumerate([(d_ray, 'Ray Theory'), (d_hffk, title_hffk)]):
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        triang = Triangulation(d['stlon'], d['stlat'])
        colors_face = cmap(norm(d['si']))
        ax.plot_trisurf(triang, d['si'], cmap=cmap, norm=norm,
                        linewidth=0, antialiased=False, alpha=0.9)
        ax.set_xlabel('Lon', fontsize=8)
        ax.set_ylabel('Lat', fontsize=8)
        ax.set_zlabel('SI', fontsize=8)
        ax.set_title(title, fontsize=11)
        ax.tick_params(labelsize=7)

    fig.suptitle('SinkingSlab 3D surface: SI(lon, lat)', fontsize=12)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=fig.axes, label='SI', shrink=0.6)
    plt.tight_layout()
    _save_or_show(fig, out)


def plot_diff_map(d_ref, d_comp, label_ref='Ray', label_comp='HFFK', out=None):
    """
    差值圖：d_comp.si - d_ref.si，顯示 HFFK 抹平的空間格局。
    """
    if len(d_ref['si']) != len(d_comp['si']):
        print("[WARN] diff: size mismatch, skipping")
        return
    diff = d_comp['si'] - d_ref['si']
    fig, ax = plt.subplots(figsize=(6, 5))
    norm = make_norm(diff)
    sc = ax.scatter(d_ref['stlon'], d_ref['stlat'], c=diff,
                    cmap=diverging_cmap(), norm=norm, s=25,
                    edgecolors='k', linewidths=0.3)
    ax.set_xlabel('Lon (°)')
    ax.set_ylabel('Lat (°)')
    ax.set_title(f'SI difference: {label_comp} − {label_ref}', fontsize=11)
    ax.set_aspect('equal', 'box')
    fig.colorbar(sc, ax=ax, label='ΔSI')
    plt.tight_layout()
    _save_or_show(fig, out)


def plot_period_std_bar(stds, periods_labels, out=None):
    """
    長條圖：各 period 的 std(SI)，顯示單調遞減趨勢。
    stds: dict mapping label → std value
    """
    fig, ax = plt.subplots(figsize=(7, 8))
    labels = list(stds.keys())
    vals = list(stds.values())
    colors = ['#555555'] + [plt.cm.Blues(0.4 + 0.12*i) for i in range(len(labels)-1)]
    bars = ax.bar(labels, vals, color=colors, edgecolor='k', linewidth=0.7)
    ax.set_ylabel('std(SI)')
    ax.set_title('Period effect: std(HFFK SI) vs period\nMonotonic decrease = Fresnel zone averaging', fontsize=10)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.0001,
                f'{v:.4f}', ha='center', va='bottom', fontsize=9)
    ax.set_ylim(0, max(vals)*1.15)
    plt.tight_layout()
    _save_or_show(fig, out)


def _save_or_show(fig, out):
    if out:
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved: {out}")
    else:
        plt.show()
    plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='SinkingSlab HFFK 3D visualization')
    parser.add_argument('--ray',       required=True, help='ray theory SI .dat file')
    parser.add_argument('--hffk3',     default=None,  help='HFFK 3s SI .dat file')
    parser.add_argument('--hffk5',     default=None,  help='HFFK 5s SI .dat file')
    parser.add_argument('--hffk8',     default=None,  help='HFFK 8s SI .dat file')
    parser.add_argument('--hffk15',    default=None,  help='HFFK 15s SI .dat file')
    parser.add_argument('--hffk25',    default=None,  help='HFFK 25s SI .dat file')
    # SplittingParameters files (for fast axis quiver overlay)
    parser.add_argument('--sp_ray',   default=None,  help='Ray theory SplittingParameters .dat')
    parser.add_argument('--sp_hffk8', default=None,  help='HFFK 8s SplittingParameters .dat')
    parser.add_argument('--sp_hffk25',default=None,  help='HFFK 25s SplittingParameters .dat')
    parser.add_argument('--quiver_skip', type=int, default=3,
                        help='Plot fast axis every N stations (default 3)')
    parser.add_argument('--quiver_fixed_len', action='store_true',
                        help='All fast axis arrows same length (direction only)')
    parser.add_argument('--receivers', default=None,
                        help='Receivers.dat path (default: auto-detect psi_input/Receivers.dat)')
    parser.add_argument('--src',   type=int, default=None,
                        help='Filter to single source ID (default: all sources)')
    parser.add_argument('--mode',  default='all',
                        choices=['map', '3d', 'surface', 'diff', 'bar', 'all'],
                        help='Which plot(s) to generate')
    parser.add_argument('--vmax',  type=float, default=None,
                        help='Symmetric colorbar limit (default: 98th percentile)')
    parser.add_argument('--out',   default=None,
                        help='Output filename prefix (e.g. "figs/sinkslab"). '
                             'Appends _map.png, _3d.png, etc. If omitted → interactive.')
    args = parser.parse_args()

    # ── load receivers ────────────────────────────────────────────────────────
    rec_path = args.receivers or find_receivers(args.ray)
    if rec_path:
        receivers = load_receivers(rec_path)
        print(f"Receivers: {len(receivers)} stations from {rec_path}")
    else:
        receivers = {}
        print("[WARN] Receivers.dat not found — spatial maps will have no coordinates")

    # ── load SI data ─────────────────────────────────────────────────────────
    d_ray = load_si(args.ray, args.src, receivers)
    n_coord = np.sum(~np.isnan(d_ray['stlon']))
    print(f"Ray theory: N={len(d_ray['si'])} obs, {n_coord} with coordinates, "
          f"mean={d_ray['si'].mean():.4f}, std={d_ray['si'].std():.4f}")

    hffk_files = {
        '3s':  args.hffk3,
        '5s':  args.hffk5,
        '8s':  args.hffk8,
        '15s': args.hffk15,
        '25s': args.hffk25,
    }
    d_hffk = {}
    for label, fpath in hffk_files.items():
        if fpath and Path(fpath).exists():
            d_hffk[label] = load_si(fpath, args.src, receivers)
            print(f"HFFK {label}: N={len(d_hffk[label]['si'])}, "
                  f"mean={d_hffk[label]['si'].mean():.4f}, "
                  f"std={d_hffk[label]['si'].std():.4f}")

    def out_path(suffix):
        if args.out is None:
            return None
        return f"{args.out}_{suffix}.png"

    # ── load SP (SplittingParameters) for quiver ─────────────────────────────
    sp_files = {'Ray Theory': args.sp_ray, 'HFFK 8s': args.sp_hffk8, 'HFFK 25s': args.sp_hffk25}
    d_sp_lookup = {}
    for label, fpath in sp_files.items():
        if fpath and Path(fpath).exists():
            d_sp_lookup[label] = load_sp(fpath, args.src, receivers)
            d = d_sp_lookup[label]
            print(f"SP {label}: N={len(d['dt'])}, dt_mean={d['dt'].mean():.3f}s, "
                  f"phi_mean={np.degrees(d['phi']).mean():.1f}°")

    # ── 2D map ───────────────────────────────────────────────────────────────
    if args.mode in ('map', 'all'):
        datasets = [d_ray] + list(d_hffk.values())
        titles   = ['Ray Theory'] + [f'HFFK {k}' for k in d_hffk]
        # match SP datasets to the same order
        sp_datasets = [d_sp_lookup.get(t) for t in titles] if d_sp_lookup else None
        plot_2d_maps(datasets, titles, sp_datasets=sp_datasets,
                     quiver_skip=args.quiver_skip,
                     quiver_fixed_len=args.quiver_fixed_len,
                     vmax=args.vmax, out=out_path('map'))

    # ── 3D scatter ───────────────────────────────────────────────────────────
    if args.mode in ('3d', 'all'):
        datasets = [d_ray] + list(d_hffk.values())
        titles   = ['Ray Theory'] + [f'HFFK {k}' for k in d_hffk]
        plot_3d_scatter(datasets, titles, vmax=args.vmax, out=out_path('3d'))

    # ── 3D surface ───────────────────────────────────────────────────────────
    if args.mode in ('surface', 'all') and d_hffk:
        # 用 最長 period 的 HFFK 對比
        last_label = list(d_hffk.keys())[-1]
        plot_3d_surface(d_ray, d_hffk[last_label],
                        title_hffk=f'HFFK {last_label}',
                        vmax=args.vmax, out=out_path('surface'))

    # ── diff map ─────────────────────────────────────────────────────────────
    if args.mode in ('diff', 'all') and d_hffk:
        for label, d in d_hffk.items():
            plot_diff_map(d_ray, d,
                          label_ref='Ray', label_comp=f'HFFK {label}',
                          out=out_path(f'diff_{label}'))

    # ── period std bar ───────────────────────────────────────────────────────
    if args.mode in ('bar', 'all'):
        stds = {'Ray Theory': d_ray['si'].std()}
        stds.update({f'HFFK {k}': v['si'].std() for k, v in d_hffk.items()})
        plot_period_std_bar(stds, list(stds.keys()), out=out_path('std_bar'))


if __name__ == '__main__':
    main()
