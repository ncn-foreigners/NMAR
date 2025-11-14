# Detect whether the RHS explicitly requests an intercept

The Formula/terms machinery already controls implicit intercepts, but we
warn users who add \`+ 1\` manually. This helper scans the parsed
language tree for an explicit scalar \`1\` so we can emit a user-facing
warning without treating \`.\` expansions or other implicit intercepts
as explicit requests.

## Usage

``` r
el_has_explicit_intercept(node)
```
