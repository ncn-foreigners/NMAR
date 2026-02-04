# Format an abridged call line for printing

Builds a concise one-line summary of the original call without
materializing large objects (e.g., full data frames). Intended for use
by print/summary methods.

## Usage

``` r
nmar_format_call_line(x)
```

## Details

Uses option \`nmar.show_call\` (default TRUE). Width can be tuned via
option \`nmar.call_width\` (default 120).
