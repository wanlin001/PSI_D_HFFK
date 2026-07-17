Sources_uniform48.dat — teleseismic SKS tier-1 catalogue
============================================================
Reference: Lin et al. (2014a,b) — Δ≈88–115°, deep focus (>50 km).

Structure: 48 = 16 azimuths × 3 tiers
  - 16 directions every 22.5° → 8 BAZ bins (45°) × 2 dirs per bin
  - Tier 0: Δ=88.0°, depth=100 km
  - Tier 1: Δ=100.0°, depth=300 km
  - Tier 2: Δ=115.0°, depth=600 km

Centroid: (121.0°E, 24.0°N) — Taiwan network reference

Job A/C read:
  ${PSI_DIR}/psi_input/Sources_uniform48.dat
  ${PSI_DIR}/psi_input/DUMMY_SI_uniform48.dat

After updating sources, rerun:
  sbatch scripts/slurm/jobA2_ray_si_uniform48.sh
  sbatch scripts/slurm/jobC_hffk_uniform48.sh

Sample (src1): Δ=88.0° BAZ=0.0° bin0
