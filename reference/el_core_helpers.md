# EL core helpers

Internal helpers for solving and post-processing the EL system.
[`el_run_solver()`](https://ncn-foreigners.ue.poznan.pl/NMAR/index.html/reference/el_run_solver.md)
orchestrates
[`nleqslv::nleqslv()`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html)
with a small, deterministic fallback ladder; `el_post_solution()`
computes masses and the point estimate with denominator guards and
optional trimming.
