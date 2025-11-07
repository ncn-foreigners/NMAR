# S3 helpers for NMAR engine objects

Lightweight, user-facing methods for engine configuration objects (class
\`nmar_engine\`). These improve discoverability and provide a consistent
print surface across engines while keeping the objects as simple lists
internally.

## Design

\- \`engine_name()\` returns a canonical identifier used across the
package (e.g., in \`nmar_result\$meta\$engine_name\`). -
\`print.nmar_engine()\` provides a concise, readable summary of the
engine configuration; engine-specific classes reuse the parent method
unless they need to override it. - \`engine_config()\` returns the
underlying configuration as a named list for programmatic inspection.
