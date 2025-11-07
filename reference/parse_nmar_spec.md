# NMAR input parsing and validation helpers

These developer-facing utilities centralize the parsing of user supplied
formulas/data and allow engines to declare small trait lists that
control validation strictness. Keeping all argument checks in one place
ensures consistent error messages and makes it easier to extend the
package with additional engines.

## Usage

``` r
parse_nmar_spec(formula, data, env = parent.frame())
```
