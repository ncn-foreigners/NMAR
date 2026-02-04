# Solver orchestration

Solver orchestration

## Usage

``` r
el_run_solver(
  equation_system_func,
  analytical_jac_func,
  init,
  final_control,
  top_args,
  solver_method,
  use_solver_jac,
  K_beta,
  K_aux,
  respondent_weights,
  N_pop,
  trace_level = 0
)
```

## Arguments

- equation_system_func:

  Function mapping parameter vector to equation residuals.

- analytical_jac_func:

  Analytic Jacobian function; may be NULL if unavailable or when forcing
  Broyden.

- init:

  Numeric vector of initial parameter values.

- final_control:

  List passed to `nleqslv::nleqslv(control = ...)`.

- top_args:

  List of top-level
  [`nleqslv::nleqslv`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html)
  args (e.g., `global`, `xscalm`).

- solver_method:

  Character; one of "auto", "newton", or "broyden".

- use_solver_jac:

  Logical; whether to pass analytic Jacobian to Newton.

- K_beta:

  Integer; number of response model parameters.

- K_aux:

  Integer; number of auxiliary constraints.

- respondent_weights:

  Numeric vector of base sampling weights.

- N_pop:

  Numeric; population total.

- trace_level:

  Integer; verbosity level.
