# Casimir-Constrained-Optimization

Mathematica 9 source code, numerical data, and figures accompanying the [Casimir optimization study on arXiv](https://arxiv.org/abs/2610.02471).

The study considers outward Casimir–Lifshitz pressure in a planar PEC–vacuum–biased GaAs–vacuum–PEC cavity. It maximizes the smaller pressure on the two conductors while accounting for uncertainty in the vacuum gaps and optical drive. This revision adds comparisons with gold reflectors, an alternative GaAs optical model, nominal and robust designs, uncertainty budgets, and a continuous-uncertainty audit. The saved results are numerical predictions for the stated models and constraints.

## Principal numerical results

At 300 K, the PEC design has equal nominal gaps of 1.35122 μm, a nominal drive parameter of η = 0.945, and an adopted GaAs thickness of 13.5 μm. For gap uncertainty of ±10 nm and drive uncertainty of ±0.005, its worst-case outward pressure is 2.25684 mPa, exceeding the corresponding far-field pressure by 0.18116 mPa. The gold-reflector comparison gives a robust pressure of 2.16428 mPa under the corresponding conditions. The full finite-slab calculation differs from the reduced model by approximately 0.00160% at the active worst-case point.

The underlying values are recorded in [`headline_results.csv`](data_mathematica9_pra/headline_results.csv), [`reflector_comparison.csv`](data_mathematica9_pra/reflector_comparison.csv), and [`finite_slab_validation.csv`](data_mathematica9_pra/finite_slab_validation.csv).

## Repository contents

| Path | Purpose |
| --- | --- |
| [`Casimir_Optimization_PRA_Mathematica9.m`](Casimir_Optimization_PRA_Mathematica9.m) | Main Mathematica 9 definitions, numerical study, validation, and exports. |
| [`Casimir_Optimization_PRA_Launcher.nb`](Casimir_Optimization_PRA_Launcher.nb) | Notebook that loads the source, runs the self-test, and starts the study. |
| [`data_mathematica9_pra/`](data_mathematica9_pra/) | Machine-readable results and the data underlying the figures. |
| [`figures_mathematica9_pra/`](figures_mathematica9_pra/) | Six current PDF figures exported by the revised script. |

The data directory includes the robust design map and pressure decomposition; quadrature and cutoff convergence checks; nominal-versus-robust and PEC-versus-gold comparisons; GaAs model sensitivity; uncertainty trade-offs, stationary candidates, and finite-slab validation. Figure CSVs are in the data directory, while the current PDFs are directly inside `figures_mathematica9_pra/`.

## Reproduce the study

Use Wolfram Mathematica 9. Keep the launcher notebook and the `.m` source in the same directory, start a fresh kernel, and evaluate the launcher's three input cells in order:

1. Load `Casimir_Optimization_PRA_Mathematica9.m`.
2. Evaluate `SelfTest[]` and confirm that it returns `True`.
3. Evaluate `RunStudy[]`.

Alternatively, from the repository root in a Mathematica session:

```wolfram
Get["Casimir_Optimization_PRA_Mathematica9.m"];
SelfTest[]
RunStudy[]
```

`RunStudy[]` generates the original optimization, the added revision analyses, and the finite-slab validation. It writes CSV files to `data_mathematica9_pra/` and the six current PDF figures to `figures_mathematica9_pra/`. The full calculation uses extensive numerical quadrature and uncertainty scans and may take substantial time. The repository includes saved outputs for inspection without rerunning it.

## Figures

The current PDFs are [slab-opacity indicators](figures_mathematica9_pra/fig2_opacity_final.pdf), [pressure curves](figures_mathematica9_pra/fig3_pressure_curves_final.pdf), [robust design map](figures_mathematica9_pra/fig4_robust_design_map_final.pdf), [pressure decomposition](figures_mathematica9_pra/fig5_decomposition_final.pdf), [optimum spectrum](figures_mathematica9_pra/fig6_spectrum_optimum_final.pdf), and [gap-asymmetry diagnostic](figures_mathematica9_pra/fig7_asymmetry_final.pdf).
