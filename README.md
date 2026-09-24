# Minimal MVFP collective-escape reproduction

This folder is a small, standalone MATLAB reproduction package for the paper's central mechanism. It generates two representative results and keeps the main parameters in one visible block so readers can modify the model and scan additional cases.

## Contents

- `run_minimal_reproduction.m`: the complete source. It contains the MVFP finite-volume solver, initial two-packet density, self-consistent feedback, potential, readout, and threshold-arrival calculation as local functions.
- `figures/minimal_core.{png,pdf}`: density evolution, self-consistent potential, feedback field, and finite-time readout.
- `figures/minimal_readout_mfpt.{png,pdf}`: final readout `P(T)` and first threshold-arrival time versus the initial left fraction `P0`.
- `results/minimal_readout_summary.csv`: compact numeric summary for the second figure.
- `results/minimal_summary.mat`: MATLAB result struct for further plotting.

The included figures are example outputs. The script recomputes them from scratch and writes fresh files into `figures/` and `results/`.

## Requirements

- MATLAB R2020a or newer is recommended.
- No project-specific MATLAB path, toolbox, input data, or external result file is required. The script uses `ode15s`, `exportgraphics`, and base MATLAB functionality.

## Run

From this directory:

```matlab
run('run_minimal_reproduction.m')
```

From PowerShell:

```powershell
matlab -batch "cd('C:/path/to/code_git'); run('run_minimal_reproduction.m')"
```

The run is deterministic because the seed is recorded in the configuration block. The deterministic PDE itself does not use random draws.

## Model and extension points

The script uses the pairwise-only feedback field

```text
F(t) = sigma1 * M1(t) + bias,
M1(t) = integral x rho(x,t) dx,
```

with effective potential

```text
U_eff(x;F) = -0.5*(u-sigma1)*x^2 + 0.25*v*x^4 - F*x.
```

The adjustable block at the top of `run_minimal_reproduction.m` controls `D`, `sigma1`, `bias`, grid resolution, observation time, `P0` lists, and the threshold. To explore other paper results, change those lists or add another call to `simulate_case` and plot the returned fields. The full research implementation with all supplementary scans remains in the parent `code/` package.

## Reproducibility scope

This repository intentionally contains only a small computational example and two representative outputs. It does not claim to regenerate every manuscript panel with the default command. The full scan code is available in the parent project package; the same model parameters and solver conventions can be used here to build additional figures.
