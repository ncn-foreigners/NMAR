# Create the NMAR delta indicator column

Adds a respondent indicator column to the data, choosing a unique
internal name (\`..nmar_delta..\` plus numeric suffixes when necessary)
so user columns are never clobbered.

## Usage

``` r
el_make_delta_column_name(data, outcome_var, respondent_mask = NULL)
```
