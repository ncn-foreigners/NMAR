# Parent S3 surface for NMAR results

Methods that apply to the parent \`nmar_result\` class and are not
specific to a particular engine (e.g., EL). Engines return a child class
(e.g., \`nmar_result_el\`) that inherits from \`nmar_result\` and may
override or extend behavior.

## Details

S3 surface for base \`nmar_result\`

Result objects expose a universal schema: - \`y_hat\`,
\`estimate_name\`, \`se\`, \`converged\`. - \`model\`: list with
\`coefficients\`, \`vcov\`, plus optional extras. - \`weights_info\`:
list with respondent weights and trimming metadata. - \`sample\`: list
with total units, respondent count, survey flag, and \`design\`. -
\`inference\`: variance metadata (\`variance_method\`, \`df\`,
diagnostic flags). - \`diagnostics\`, \`meta\`, and \`extra\` for
estimator-specific details.

New engines should populate these components in their constructors and
rely on the \`nmar_result_get\_\*\` utilities when implementing
child-specific S3 methods.
