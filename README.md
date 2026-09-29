# Casimir-Constrained-Optimization
# Robust optimization of nonequilibrium Casimir repulsion

Mathematica 9 source, numerical data, and figures for a constrained optimization of outward Casimir–Lifshitz pressure in a planar **PEC–vacuum–biased GaAs–vacuum–PEC** cavity. The objective is to maximize the smaller outward pressure on the two conductors under gap and drive uncertainty. In the optically thick slab model, identical boundaries admit a symmetric optimizer, reducing the robust search to one nominal gap.

## Main result

At **300 K**, the optimized nominal gaps are **1.35122 μm** on both sides, the nominal drive parameter is **η = 0.945**, and the adopted GaAs slab thickness is **13.5 μm**. With independent gap deviations of **±10 nm** and a drive deviation of **±0.005**, the worst-case outward pressure is **2.25684 mPa**. This exceeds the corresponding far-field pressure by **0.18116 mPa**. A separate full finite-slab scattering calculation differs from the reduced model by about **0.00160%** at the active worst-case point.

These are numerical results for the implemented model, material parameters, and uncertainty bounds. See [headline results](data_mathematica9_final/headline_results.csv) and [finite-slab validation](data_mathematica9_final/finite_slab_validation.csv).

## Contents

| Path | Contents |
| --- | --- |
| [`Casimir_Optimization_Mathematica9_Final_2.m`](Casimir_Optimization_Mathematica9_Final_2.m) | Main definitions, optimization, checks, and data/figure exports. |
| [`Casimir_Optimization_Mathematica9_Final_2.nb`](Casimir_Optimization_Mathematica9_Final_2.nb) | Notebook version of the main calculation. |
| [`data_mathematica9_final/`](data_mathematica9_final/) | Saved results, figure data, and convergence/uncertainty checks. |
| [`figures_mathematica9_final/Clean_Figs/`](figures_mathematica9_final/Clean_Figs/) | Final figure PDFs, CSVs, and figure-preparation notebook. |

The calculation combines a GaAs optical response, equilibrium Lifshitz pressure, and propagating and evanescent nonequilibrium contributions. It scans the design range, refines the best feasible gap, validates the uncertainty minimum on successively finer grids, and exports the results.

## Reproduce

In **Wolfram Mathematica 9**, start in the repository root and evaluate:

```wolfram
Get["Casimir_Optimization_Mathematica9_Final_2.m"];
SelfTest[]
RunStudy[]
RunFiniteSlabValidation[]
```

Check that `SelfTest[]` returns `True` before running the full study. `RunStudy[]` writes CSVs to `data_mathematica9_final/` and computational figure PDFs to `figures_mathematica9_final/`. `RunFiniteSlabValidation[]` separately writes `finite_slab_validation.csv`. The full numerical study can take substantial time.

For the final styled figure PDFs in `figures_mathematica9_final/Clean_Figs/`, use the [figure-preparation notebook](figures_mathematica9_final/Clean_Figs/Clean_Final_Figs_Based_on_Casimir_Optimization_Mathematica9_6.nb). The main script exports computational figure PDFs to its own output directory; the styled versions are separate files in `Clean_Figs/`.
