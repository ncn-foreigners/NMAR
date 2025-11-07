# Create Verbose Printer Factory

Creates a verbose printing function based on trace level settings.
Messages are printed only if their level is \<= trace_level.

## Usage

``` r
create_verboser(trace_level = 0)
```

## Arguments

- trace_level:

  Integer 0-3; controls verbosity detail: - 0: No output (silent mode) -
  1: Major steps only (initialization, convergence) - 2: Moderate detail
  (iteration summaries, key diagnostics) - 3: Full detail (all
  diagnostics, intermediate values)

## Value

A function with signature: \`verboser(msg, level = 1, type = c("info",
"step", "detail", "result"))\`
