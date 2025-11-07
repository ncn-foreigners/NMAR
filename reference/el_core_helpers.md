# EL core helpers

Internal helpers for solving and post-processing the EL system.
\`el_run_solver()\` orchestrates \`nleqslv\` with a small, deterministic
fallback ladder; \`el_post_solution()\` computes masses and the point
estimate with denominator guards and optional trimming.
