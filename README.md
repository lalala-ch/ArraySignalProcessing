# ArraySignalProcessing

阵列信号处理的 MATLAB 学习笔记与仿真代码，包含窄带波束形成和窄带波达方向（DOA）估计。Notebook 将原理推导、仿真代码与结果放在一起，独立 `.m` 文件提供求解函数或补充示例。

## 目录

```text
ArraySignalProcessing/
├── narrowband beamforming/       # 波束形成：7 个 Notebook、5 个示例脚本
├── narrowband DOA estimation/    # DOA 估计：2 个 Notebook、4 个求解函数
├── README.md
└── LICENSE
```

## Narrowband beamforming

| Notebook | 包含的算法或内容 |
| --- | --- |
| [Conventional_Beamformer.ipynb](narrowband%20beamforming/Conventional_Beamformer.ipynb) | 常规波束形成、阵列方向图 |
| [FIR_beampattern.ipynb](narrowband%20beamforming/FIR_beampattern.ipynb) | 利用 FIR 滤波器设计方法构造阵列方向图 |
| [LSbasedBeamformer.ipynb](narrowband%20beamforming/LSbasedBeamformer.ipynb) | 最小二乘波束形成与方向图拟合 |
| [MVDR_LCMVBeamformer.ipynb](narrowband%20beamforming/MVDR_LCMVBeamformer.ipynb) | MVDR / Capon、线性约束最小方差（LCMV）波束形成 |
| [RefsigBasedBeamformer.ipynb](narrowband%20beamforming/RefsigBasedBeamformer.ipynb) | 基于参考信号的 LMS、RLS 自适应波束形成 |
| [BlindAdaptiveBeamformer.ipynb](narrowband%20beamforming/BlindAdaptiveBeamformer.ipynb) | 最小功率法、CMA、LS-CMA、LS-SCORE 盲波束形成 |
| [RobustBF.ipynb](narrowband%20beamforming/RobustBF.ipynb) | 对角加载 Capon、最坏情况鲁棒波束形成、Gu–Leshem 干扰协方差重构与导向向量估计 |

补充脚本：

- `LSbasedBeamformer_2.m`：最小二乘波束形成示例。
- `qpsk_ls_score_demo_fixed.m`：QPSK 信号的 LS-SCORE 示例。
- `wc_rab_demo.m`：利用等效对角加载与二分搜索求解最坏情况鲁棒波束形成，不依赖 CVX。
- `wc_rab_cvx_demo.m`：使用 CVX 求解最坏情况鲁棒波束形成。
- `gu_leshem_capon_reconstruction_demo.m`：Gu–Leshem 协方差重构与导向向量估计示例，依赖 CVX。

## Narrowband DOA estimation

| Notebook | 包含的算法或内容 |
| --- | --- |
| [Classical_NB_Doa.ipynb](narrowband%20DOA%20estimation/Classical_NB_Doa.ipynb) | Bartlett、Capon、FFT 空间谱估计；MUSIC、FFT 加速 MUSIC、LS-ESPRIT、TLS-ESPRIT；确定性最大似然（DML）的交替投影、Gauss–Newton 和 EM 求解 |
| [SparsityBased_NB_Doa.ipynb](narrowband%20DOA%20estimation/SparsityBased_NB_Doa.ipynb) | OMP、CoSaMP；组稀疏优化的 FISTA、ADMM；ℓ1-SVD；FOCUSS；SPICE 稀疏协方差拟合 |

独立求解函数与 Notebook 位于同一目录：

- `solve_group_lasso_fista.m`：FISTA 组稀疏求解，支持自适应重启。
- `solve_group_lasso_admm.m`：ADMM 组稀疏求解，记录原始/对偶残差和停止阈值。
- `solve_focuss.m`：FOCUSS 迭代求解。
- `solve_spice.m`：统一处理 SCM 满秩与秩亏情况，估计角度网格功率和各阵元噪声功率，并返回拟合协方差及迭代记录。

## 运行方法

1. 安装 MATLAB，以及支持 MATLAB 的 Jupyter 内核。本仓库 Notebook 的内核名为 `jupyter_matlab_kernel`，不是 Python 内核；部分求解函数已在 MATLAB R2024b 下验证。
2. 按示例需要安装工具箱：Signal Processing Toolbox（如 `findpeaks`、`firpm`、`fir1`、`upfirdn`）和 Communications Toolbox（如 `awgn`、`pskmod`、`rcosdesign`）。使用 CVX 的鲁棒波束形成示例还需单独安装 CVX，并先运行 `cvx_setup`。
3. 在对应子目录中打开 Notebook，将 MATLAB 工作目录切换到该子目录，按单元格顺序运行；独立 `.m` 示例可以直接在 MATLAB 中运行。
4. 稀疏 DOA Notebook 需要同目录下的四个求解函数。也可以从仓库根目录设置路径：

   ```matlab
   addpath(fullfile(pwd, 'narrowband beamforming'));
   addpath(fullfile(pwd, 'narrowband DOA estimation'));
   ```

示例使用模拟数据，部分算法共用前面的初始化结果，请勿跳过数据生成单元格。随机数据、参数和迭代上限会影响结果；达到迭代上限不代表满足收敛阈值。仓库保留原 Notebook 中已有的输出，但不表示所有示例都已在当前环境重新运行验证。

## 许可证

本项目采用 [MIT License](LICENSE)。算法及论文的相关说明和引用保留在 Notebook 与示例脚本中。
