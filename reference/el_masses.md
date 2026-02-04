# Compute probability masses

Compute probability masses

## Usage

``` r
el_masses(weights, denom, floor, trim_cap)
```

## Arguments

- weights:

  numeric respondent base weights (d_i)

- denom:

  numeric denominators Di after floor guard

- floor:

  numeric small positive guard

- trim_cap:

  numeric cap (\>0) or Inf

## Value

list with mass_untrim, mass_trimmed, prob_mass, trimmed_fraction
