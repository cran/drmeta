# Cross-package validation

`drmeta_numerical_validation.R` records the independent numerical audit used
to validate the package. The likelihood audit was initially run with R 4.3.3,
metafor 4.4-0, metadat 1.2-0, and drmeta 0.2.2. The complete script was rerun
successfully through the drmeta 0.2.3 finalization gate on 2026-09-28, and the
principal numerical invariants are also conditional testthat tests.

The script checks ML estimates and likelihoods and verifies REML point
estimates while retaining the implementations' different likelihood constants,
checks the nested joint model, documents bounded-optimizer start sensitivity,
fits the categorical scale diagnostic, and reproduces the influence and score
sensitivity summaries used in the article.
