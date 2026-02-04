# Run method for EL engine

Run method for EL engine

## Usage

``` r
# S3 method for class 'nmar_engine_el'
run_engine(engine, formula, data, trace_level = 0)
```

## Arguments

- engine:

  An object of class `nmar_engine_el`.

- formula:

  A two-sided formula passed through by
  [`nmar()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/nmar.md).

- data:

  A `data.frame` or `survey.design`.

- trace_level:

  Integer 0-3 controlling verbosity.

## Value

An object of class `nmar_result_el`.
