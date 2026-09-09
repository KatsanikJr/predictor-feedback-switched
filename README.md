# Predictor-Feedback Stabilization of Linear Switched Systems with State-Dependent Switching or Delayed Switching Input

The purpose of these MATLAB codes is to reproduce the numerical examples presented in Section V of the paper

**"Predictor-Feedback Stabilization of Linear Switched Systems with State-Dependent Switching or Delayed Switching Input"**

by Andreas Katsanikakis, Nikolaos Bekiaris-Liberis, and Delphine Bresch-Pietri.

The paper considers predictor-feedback stabilization for two classes of switched linear systems: systems with state-dependent switching and a delayed continuous control input, and systems in which the switching input itself is delayed.

## Requirements

The simulation codes require:

- **MATLAB**
- **CVX toolbox** only for running `find_LMI_IV_B.m`

The main simulation scripts do not require CVX.

## Usage

The repository contains the MATLAB scripts corresponding to the two numerical examples in Section V of the paper.

### Section V.A — State-dependent switching with delayed control input

This example illustrates Theorem 3.1.

Run:

```matlab
Sim_IV_A.m
```

The simulation data are saved in `myworkspace_IV_A.mat`.

Then run:

```matlab
plots_IV_A.m
```

to generate the corresponding state and control-input plots.

### Section V.B — Delayed switching input

This example illustrates Theorem 4.3.

Run:

```matlab
Sim_IV_B.m
```

The simulation data are saved in `myworkspace_IV_B.mat`.

Then run:

```matlab
plots_IV_B.m
```

to generate the corresponding state and switching-input plots.

The file

```matlab
find_LMI_IV_B.m
```

contains the grid-based BMI/LMI search used to obtain the feasible Lyapunov matrices employed in the numerical example of Section V.B. This script requires CVX and is not needed to run the simulation itself.

## Files

- `Sim_IV_A.m` — simulation of the example in Section V.A
- `plots_IV_A.m` — plots for the example in Section V.A
- `Sim_IV_B.m` — simulation of the example in Section V.B
- `plots_IV_B.m` — plots for the example in Section V.B
- `find_LMI_IV_B.m` — auxiliary BMI/LMI search for the matrices used in Section V.B

## License

Copyright Andreas Katsanikakis 2026. See LICENSE.txt for licensing information.

## Acknowledgements

Funded by the European Union (ERC, C-NORA, 101088147). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Cite this work

If you use these codes, please cite:

@unpublished{katsbekbrp,
  title = {Predictor-Feedback Stabilization of Linear Switched Systems with State-Dependent Switching or Delayed Switching Input},
  author = {Katsanikakis, Andreas and Bekiaris-Liberis, Nikolaos and Bresch-Pietri, Delphine},
  url = {https://hal.science/hal-05579119},
  note = {preprint},
  year = {2026},
  month = September,
}

Full bibliographic information will be added once available.
