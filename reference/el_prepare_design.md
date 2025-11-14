# EL design construction and Formula workflow

Parses the user formula once via \`model.frame()\` on a \`Formula\`
object and materializes each RHS block with the standard
\`model.matrix()\` machinery. Auxiliary designs are built on the full
dataset (with the intercept and outcome columns removed) while the
missingness design uses only respondent rows and always includes an
explicit intercept and outcome column.

## Usage

``` r
el_prepare_design(formula, data)
```
