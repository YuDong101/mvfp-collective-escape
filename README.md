# MVFP collective escape

Numerical and plotting code for *From slow redistribution to collective escape: self-consistent feedback and state selection*, by Dong Yu, Hang Chen, Qiming Pei, Long Jiang and Ya Jia. Correspondence: yudong@yangtzeu.edu.cn; jiay@ccnu.edu.cn.

The package reproduces Figures 1–5 and the supporting local theory. Figure 6 calculations can be obtained by adapting the supplied density and particle routines to the parameters below. Its numerical data are not included.

## Requirements

The calculations were tested with MATLAB R2024b. The default numerical stages run serially. Plotting was tested with Python 3.13.13; the package versions are pinned in `requirements-plotting.txt`.

```powershell
python -m pip install -r requirements-plotting.txt
```

On Windows, extract the code to a short path such as `C:\MVFP_reproduction`. The inherited checkpoint filenames contain full hashes, and deeply nested paths can exceed MATLAB/HDF5 filename limits.

## Quick numerical check

From the code directory in MATLAB:

```matlab
run('verify_reproduction.m')
```

This performs new small calculations of upwind and fitted-flux density evolution, frozen-feedback controls, finite-particle occupancy, stationary branches, free energy, time-window readouts and an unstable eigenmode. It also checks MATLAB dependencies. Results are written to `verification/quick/`. No precomputed study data are needed for this check.

## Full numerical reproduction

```matlab
run('run_recomputation.m')
```

The switches at the top select stationary branches/free-energy trajectories, Case A, Case B and final data assembly. All are enabled by default. The production settings are retained; the full run takes substantially longer than the quick check.

Each stage writes its own outputs. Disabled stages retain existing files. Compatible checkpoints let interrupted calculations resume. Use a fresh code/output directory after changing the model, numerical settings or source code; do not combine checkpoints from the earlier two-directory package with this layout.

After all numerical stages finish:

```powershell
python recompute/analyze_boundary.py
python render_fresh.py --figures 1,2,3,4,5
python recompute/audit_fresh.py
```

Figure arrays are assembled under `data/fig123/` and `data/fig45/`. Internal theoretical outputs are under `outputs/pairwise_paper_revision_v2/`. Rendered PDF, PNG and SVG files go to `figures/`; the numerical audit goes to `provenance/final_numerical_audit.json`.

## Rendering the supplied Source Data

The separate `MVFP_Source_Data_fresh_2026-10-10.zip` contains the newly calculated Figure 1–5 data and supporting theoretical records from 10 October 2026. Extract it, copy its `data/` directory into this code directory, then run:

```powershell
python render_fresh.py --figures 1,2,3,4,5
```

The data package's `SOURCE_DATA_GUIDE.md` lists the scientific arrays and variable definitions. Rendering from Source Data does not rerun the simulations. The full numerical audit additionally expects the stage logs and internal output layout created by `run_recomputation.m`.

## Model and numerical methods

All model quantities and physical times are dimensionless. Case A uses `(D,sigma1)=(0.005,0.5)` and Case B uses `(0.01,2)`. Both use `u=4`, `v=20`, `bias=-0.5`, `sigma2=0`, density domain `[-1,1]` and initial packet width `0.02`. Occupancy is the density integral over `x<=0`, or the fraction of particles with `x<=0`. The instantaneous saddle is used for potential-barrier diagnostics.

`function/` is the single shared MATLAB library. Its density integrator uses the original `ode15s` branch for upwind fluxes and the validated sparse bordered formulation for pairwise exponential-fit fluxes. Fixed-feedback fitted-flux controls use the corresponding tridiagonal Jacobian. The broader free-energy illustrations retain their upwind settings; controls and local theory retain their fitted-flux settings, grids and tolerances.

The source data include comparisons with newly computed reference trajectories on 2,001- and 4,001-point grids. Numerical refinement and prediction intervals are empirical allowances, not statistical confidence intervals or rigorous continuum error bounds. `wallSeconds`, when present, measures computer execution time.

## Figure 6 parameter adaptations

Use `solve_mvfp_density_trace` for density evolution and `langevin_final_left_fraction` or `langevin_trace_realization` for finite-particle realizations. For example, a single density/particle endpoint comparison can be computed as follows:

```matlab
addpath('function');
base=mvfp_paper_revision_config(pwd);
par=base.parA; P0=0; T=100;
solver=base.mvfp;
solver.fluxScheme='exponential-fit';
solver.xGridN=1001;
trace=solve_mvfp_density_trace(par,P0,0:0.2:T,solver);
particle=struct('N',1000,'nRealizations',30,'dt',0.02,'xClip',8);
[Pmean,Pstd,values]=langevin_final_left_fraction(par,P0,T,particle,11);
stderr=Pstd/sqrt(particle.nRealizations);
```

Adapt the model parameters and initial weights to the following protocol:

- Time-series panels: `sigma1=0.5`, `D=[0.001,0.005,0.01,0.05]`, `P0=[0,0.2,0.4,0.8]`, with recorded times `0:0.2:100`.
- Response maps: use the exact `fig6_D` and `fig6_sigma1` arrays in `provenance/input_config.json`, with `P0=[0,0.1,0.2,0.3,0.4,0.5]`. The particle endpoint comparison uses `P0=0`.
- Density settings: 1,001 points on `[-1,1]`, exponential-fit fluxes, nodal control volumes, fixed readout threshold, relative tolerance `1e-6`, absolute tolerance `1e-10`. These settings are inherited from `base.mvfp` with the flux choice above. Record raw mass, boundary mass and minimum density.
- Particle settings: 1,000 particles, 30 realizations, Euler–Maruyama step `0.02` and clipping bound `8`. Initial particle packets are untruncated Gaussians; the density grid has finite endpoints.
- To reproduce the scan's random-stream convention, pass endpoint seed `11+1000*iD+iSigma`; realizations then advance the same stream sequentially. For time-series realizations, use `11+50000*iD+1000*iP+r`. Indices start at 1 and `r=1,...,30`.

Save individual realization occupancies and seeds before calculating means, sample standard deviations and standard errors. These instructions do not represent a completed new Figure 6 dataset.

## Files and provenance

| Location | Purpose |
|---|---|
| `run_recomputation.m` | Full numerical entry |
| `verify_reproduction.m`, `tests/` | Independent small numerical checks |
| `function/` | Shared solver and analysis functions |
| `recompute/` | Stage drivers, boundary prediction and audits |
| `render_fresh.py`, `plotting/` | Figure rendering |
| `provenance/input_config.json` | Model parameters and Figure 6 scan settings |
| `provenance/source_manifest.json` | Hashes of the distributed code and documentation |
| `provenance/origin_map.json` | Mapping to the executed source snapshots |

The source-integrity check can be run independently:

```powershell
python recompute/audit_sources.py
```

The compact layout shares identical functions and excludes unused higher-order-interaction modules, historical MATLAB plot helpers and duplicated wrappers. Calculated data, logs and plots are excluded from the code archive.

License selection remains pending the authors' decision. No public release or DOI is asserted by this local package.
