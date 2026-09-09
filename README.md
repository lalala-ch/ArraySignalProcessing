# ArraySignalProcessing

MATLAB study notes and simulations for array signal processing, covering narrowband beamforming and direction-of-arrival (DOA) estimation. The notebooks combine theoretical derivations, simulation code, and results. Only the `.m` solver functions required by the notebooks are included; standalone example scripts and `.mlx` files are excluded.

## Repository structure

```text
ArraySignalProcessing/
├── narrowband beamforming/       # 7 notebooks
├── narrowband DOA estimation/    # 2 notebooks and 4 solver functions
├── README.md
└── LICENSE
```

## Narrowband beamforming

| Notebook | Algorithms and topics |
| --- | --- |
| [Conventional_Beamformer.ipynb](narrowband%20beamforming/Conventional_Beamformer.ipynb) | Conventional beamforming and array beampatterns |
| [FIR_beampattern.ipynb](narrowband%20beamforming/FIR_beampattern.ipynb) | Array beampattern synthesis using FIR filter design methods |
| [LSbasedBeamformer.ipynb](narrowband%20beamforming/LSbasedBeamformer.ipynb) | Least-squares beamforming and beampattern fitting |
| [MVDR_LCMVBeamformer.ipynb](narrowband%20beamforming/MVDR_LCMVBeamformer.ipynb) | MVDR / Capon and linearly constrained minimum variance (LCMV) beamforming |
| [RefsigBasedBeamformer.ipynb](narrowband%20beamforming/RefsigBasedBeamformer.ipynb) | Reference-signal-based LMS and RLS adaptive beamforming |
| [BlindAdaptiveBeamformer.ipynb](narrowband%20beamforming/BlindAdaptiveBeamformer.ipynb) | Minimum-power, CMA, LS-CMA, and LS-SCORE blind beamforming |
| [RobustBF.ipynb](narrowband%20beamforming/RobustBF.ipynb) | Diagonally loaded Capon, worst-case robust beamforming, and Gu–Leshem interference covariance reconstruction with steering vector estimation |

## Narrowband DOA estimation

| Notebook | Algorithms and topics |
| --- | --- |
| [Classical_NB_Doa.ipynb](narrowband%20DOA%20estimation/Classical_NB_Doa.ipynb) | Bartlett, Capon, and FFT spatial spectrum estimation; MUSIC, FFT-accelerated MUSIC, LS-ESPRIT, and TLS-ESPRIT; deterministic maximum likelihood (DML) estimation using alternating projection, Gauss–Newton, and EM |
| [SparsityBased_NB_Doa.ipynb](narrowband%20DOA%20estimation/SparsityBased_NB_Doa.ipynb) | OMP and CoSaMP; FISTA and ADMM for group-sparse optimization; ℓ1-SVD; FOCUSS; SPICE sparse covariance fitting |

Standalone solver functions are located alongside the notebooks:

- `solve_group_lasso_fista.m`: A FISTA group-sparse solver with adaptive restart.
- `solve_group_lasso_admm.m`: An ADMM group-sparse solver that records primal/dual residuals and stopping thresholds.
- `solve_focuss.m`: An iterative FOCUSS solver.
- `solve_spice.m`: A unified solver for full-rank and rank-deficient sample covariance matrices (SCMs). It estimates angular-grid powers and per-sensor noise powers, and returns the fitted covariance matrix and iteration histories.

## Getting started

1. Install MATLAB and a MATLAB-compatible Jupyter kernel. These notebooks use the kernel name `jupyter_matlab_kernel`, not a Python kernel. Selected solver functions have been verified in MATLAB R2024b.
2. Install the toolboxes needed by each example: Signal Processing Toolbox (e.g., `findpeaks`, `firpm`, `fir1`, and `upfirdn`) and Communications Toolbox (e.g., `awgn`, `pskmod`, and `rcosdesign`). Robust beamforming examples that use CVX require a separate CVX installation; run `cvx_setup` before using them.
3. Open a notebook in its subdirectory, set the MATLAB working directory to that subdirectory, and run the cells in order.
4. The sparsity-based DOA notebook requires the four solver functions in the same directory. Alternatively, add both subdirectories to the MATLAB path from the repository root:

   ```matlab
   addpath(fullfile(pwd, 'narrowband beamforming'));
   addpath(fullfile(pwd, 'narrowband DOA estimation'));
   ```

The examples use simulated data, and some algorithms share initialization results from earlier cells. Do not skip the data-generation cells. Results depend on random data, parameter choices, and iteration limits; reaching the iteration limit does not imply that the convergence threshold has been met. Existing notebook outputs are preserved, but not every example has been rerun and verified in the current environment.

## License

This project is released under the [MIT License](LICENSE). Algorithm descriptions and references to relevant papers are retained in the notebooks.
