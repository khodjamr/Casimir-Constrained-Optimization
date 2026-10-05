# Casimir-Constrained-Optimization

Machine-readable data and analysis files accompanying the study [arXiv:2610.02471](http://arxiv.org/abs/2610.02471).

This repository provides all the numerical data underlying Figs.2-7 and the reported validation
tables, together with the Wolfram Mathematica source code used to generate them.

## Overview

This work shows how to design a light-driven microscopic cavity that maintains a repulsive Casimir force despite realistic fabrication and operating uncertainties.

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
