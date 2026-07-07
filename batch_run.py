#!/usr/bin/env python3
"""
batch_run.py — PSI_D / ECOMAN 批次出圖
=======================================
用法：
  python batch_run.py [--config batch_config.yaml] [--model SinkingSlab] [--dry-run]

每個 model 產出：
  figs/<model>_map.png        — SI 2D map + phi tick mark + 1σ fan
  figs/<model>_misfit.png     — phi misfit map (syn − obs)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import LinearNDInterpolator
import yaml

# ── 引入 plot_si_3d 的 loader 函數 ──────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "validation"))
from plot_si_3d import (
    load_receivers, find_receivers,
    load_si, load_sp,
    plot_2d_maps,
    diverging_cmap, make_norm,
    _save_or_show,
)


# ============================================================
# Observation loader
# ============================================================
def load_obs_sp(path, receivers=None):
    """讀 PSI_D SplittingParameters 格式的 observations（同 load_sp 但不做 circular mean）。
    回傳 list of dict: {station, lon, lat, phi_deg, dt, std_phi_deg, std_dt}
    """
    records = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = [x.strip() for x in ln.split(',')]
            try:
                dt      = float(parts[0])
                phi_rad = float(parts[1])
                dt_err  = float(parts[2])
                phi_err = float(parts[3])
                station = parts[7]
            except (ValueError, IndexError):
                continue
            lon, lat = np.nan, np.nan
            if receivers and station in receivers:
                lon, lat = receivers[station]
            records.append(dict(
                station=station, lon=lon, lat=lat,
                phi_deg=np.degrees(phi_rad), dt=dt,
                std_phi_deg=np.degrees(phi_err), std_dt=dt_err,
            ))
    return [r for r in records if np.isfinite(r['lon'])]


# ============================================================
# Misfit helpers (adapted from plot_slab_contour_drift.py)
# ============================================================
def _signed_axial_diff_deg(model_phi, obs_phi):
    return (np.asarray(model_phi) - np.asarray(obs_phi) + 90.0) % 180.0 - 90.0


def _model_uncertainty_min_misfit_deg(model_phi, obs_phi, model_std_phi):
    raw = float(_signed_axial_diff_deg(model_phi, obs_phi))
    if not np.isfinite(model_std_phi):
        return raw
    tol = max(0.0, float(model_std_phi))
    return float(np.sign(raw) * max(0.0, abs(raw) - tol))


def _build_model_interpolators(sp):
    """sp = dict from load_sp: stlon, stlat, phi(rad), phi_std(rad), dt"""
    valid = ~np.isnan(sp['stlon'])
    if valid.sum() < 3:
        return None
    pts = np.column_stack([sp['stlon'][valid], sp['stlat'][valid]])
    phi_rad = sp['phi'][valid]
    phi2 = 2 * phi_rad
    return {
        'cos2phi': LinearNDInterpolator(pts, np.cos(phi2), fill_value=np.nan),
        'sin2phi': LinearNDInterpolator(pts, np.sin(phi2), fill_value=np.nan),
        'dt':      LinearNDInterpolator(pts, sp['dt'][valid], fill_value=np.nan),
        'std_phi': LinearNDInterpolator(pts, np.degrees(sp['phi_std'][valid]), fill_value=np.nan),
    }


def plot_misfit_map(sp_datasets, titles, obs_records,
                   misfit_mode='model_phi_uncertainty_min',
                   misfit_vmin=-45, misfit_vmax=45, out=None):
    """
    3 panels (one per method): scatter circle colored by phi misfit at obs stations.
    obs bar overlay shows observed phi/dt.
    """
    n = len(sp_datasets)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4.8))
    if n == 1:
        axes = [axes]

    cmap = plt.get_cmap('RdBu_r')
    norm = mcolors.Normalize(vmin=misfit_vmin, vmax=misfit_vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    for ax, sp, title in zip(axes, sp_datasets, titles):
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Lon (°)')
        ax.set_ylabel('Lat (°)')
        ax.set_aspect('equal', 'box')

        if sp is None or not obs_records:
            ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center')
            continue

        interp = _build_model_interpolators(sp)
        if interp is None:
            ax.text(0.5, 0.5, 'interp failed', transform=ax.transAxes, ha='center')
            continue

        n_drawn = 0
        for obs in obs_records:
            pt = np.array([[obs['lon'], obs['lat']]])
            cos2 = float(interp['cos2phi'](pt)[0])
            sin2 = float(interp['sin2phi'](pt)[0])
            std  = float(interp['std_phi'](pt)[0])
            if not np.isfinite(cos2) or not np.isfinite(sin2):
                continue
            model_phi = (0.5 * np.degrees(np.arctan2(sin2, cos2))) % 180.0
            if misfit_mode == 'model_phi_uncertainty_min':
                dphi = _model_uncertainty_min_misfit_deg(model_phi, obs['phi_deg'], std)
            else:
                dphi = float(_signed_axial_diff_deg(model_phi, obs['phi_deg']))

            # observed bar
            half = obs['dt'] * 0.35 / 2.0
            az = np.radians(obs['phi_deg'])
            cos_lat = np.cos(np.radians(obs['lat']))
            dlon = np.sin(az) / cos_lat * half
            dlat = np.cos(az) * half
            ax.plot([obs['lon'] - dlon, obs['lon'] + dlon],
                    [obs['lat'] - dlat, obs['lat'] + dlat],
                    'k-', lw=0.8, alpha=0.7, solid_capstyle='round')
            # misfit circle
            ax.scatter([obs['lon']], [obs['lat']],
                       s=35, marker='o', facecolors=[cmap(norm(dphi))],
                       edgecolors='k', linewidths=0.3, alpha=0.95, zorder=5)
            n_drawn += 1

        ax.text(0.02, 0.98, f'N_obs={n_drawn}',
                transform=ax.transAxes, va='top', fontsize=8)

    fig.suptitle('phi misfit: syn − obs (model uncertainty min)', fontsize=12)
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.025])
    fig.colorbar(sm, cax=cbar_ax, orientation='horizontal',
                 label='phi misfit (deg)  [model incl. ±std − obs]')
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    _save_or_show(fig, out)


# ============================================================
# Config loader
# ============================================================
def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_path(base_dir, p):
    if p is None:
        return None
    p = Path(p)
    return p if p.is_absolute() else base_dir / p


# ============================================================
# Main
# ============================================================
def run_model(model_cfg, global_cfg, base_dir, dry_run=False):
    name = model_cfg['name']
    print(f"\n{'='*60}")
    print(f"Model: {name} — {model_cfg.get('label', '')}")
    print(f"{'='*60}")

    # merged plot settings
    plot_cfg = {**global_cfg.get('plot', {}), **model_cfg.get('plot', {})}
    out_dir = base_dir / plot_cfg.get('out_dir', 'figs')
    out_dir.mkdir(parents=True, exist_ok=True)

    # receivers
    rec_path = resolve_path(base_dir, 'psi_input/Receivers.dat')
    receivers = load_receivers(str(rec_path)) if rec_path and rec_path.exists() else {}
    print(f"  Receivers: {len(receivers)} stations")

    if dry_run:
        print("  [dry-run] skipping actual plot")
        return

    # load observations
    obs_records = []
    for obs_cfg in global_cfg.get('observations', []):
        obs_path = resolve_path(base_dir, obs_cfg['path'])
        if obs_path and obs_path.exists():
            recs = load_obs_sp(str(obs_path), receivers)
            print(f"  Obs [{obs_cfg['label']}]: {len(recs)} records")
            obs_records.extend(recs)
        else:
            print(f"  Obs [{obs_cfg.get('label','')}]: not found — {obs_path}")

    # load SI and SP for each method
    si_datasets, sp_datasets, titles = [], [], []
    for m in model_cfg.get('methods', []):
        label = m['label']
        mtype = m['type']

        if mtype == 'ecoman':
            print(f"  [{label}] ECOMAN SKS-SPLIT: not yet implemented in this plotter")
            continue

        si_path = resolve_path(base_dir, m.get('si_file'))
        sp_path = resolve_path(base_dir, m.get('sp_file'))

        if si_path and si_path.exists():
            d_si = load_si(str(si_path), src=None, receivers=receivers)
            print(f"  [{label}] SI: N={len(d_si['si'])}, std={d_si['si'].std():.4f}")
        else:
            print(f"  [{label}] SI: not found — {si_path}")
            d_si = None

        if sp_path and sp_path.exists():
            d_sp = load_sp(str(sp_path), src=None, receivers=receivers)
            print(f"  [{label}] SP: N={len(d_sp['dt'])}, dt_mean={d_sp['dt'].mean():.3f}s")
        else:
            print(f"  [{label}] SP: not found — {sp_path}")
            d_sp = None

        if d_si is not None:
            si_datasets.append(d_si)
            sp_datasets.append(d_sp)
            titles.append(label)

    if not si_datasets:
        print("  No data loaded, skipping")
        return

    # ── 2D map ───────────────────────────────────────────────
    map_out = str(out_dir / f"{name}_map.png")
    plot_2d_maps(si_datasets, titles,
                 sp_datasets=sp_datasets,
                 quiver_skip=plot_cfg.get('quiver_skip', 1),
                 quiver_fixed_len=plot_cfg.get('quiver_fixed_len', True),
                 vmax=plot_cfg.get('vmax_si'),
                 out=map_out)

    # ── misfit map ───────────────────────────────────────────
    if obs_records and any(s is not None for s in sp_datasets):
        misfit_out = str(out_dir / f"{name}_misfit.png")
        plot_misfit_map(
            sp_datasets, titles, obs_records,
            misfit_mode=plot_cfg.get('misfit_mode', 'model_phi_uncertainty_min'),
            misfit_vmin=plot_cfg.get('misfit_vmin', -45),
            misfit_vmax=plot_cfg.get('misfit_vmax',  45),
            out=misfit_out,
        )
    else:
        print("  Skipping misfit map (no obs or no SP data)")


def main():
    parser = argparse.ArgumentParser(description='PSI_D batch plot runner')
    parser.add_argument('--config', default='batch_config.yaml')
    parser.add_argument('--model',  default=None, help='只跑指定 model name')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path(__file__).parent / config_path
    base_dir = config_path.parent

    cfg = load_config(config_path)
    global_cfg = cfg.get('global', {})
    models = cfg.get('models', [])

    if args.model:
        models = [m for m in models if m['name'] == args.model]
        if not models:
            print(f"Model '{args.model}' not found in config")
            sys.exit(1)

    for model_cfg in models:
        run_model(model_cfg, global_cfg, base_dir, dry_run=args.dry_run)

    print("\nDone.")


if __name__ == '__main__':
    main()
