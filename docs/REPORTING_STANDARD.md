# Equation-linked reproduction reporting standard

Each report follows **source equation → code → conditions → measured result → interpretation**.

1. Identify paper title/authors/year/DOI, exact equation-numbering version, pages and figure scope. Distinguish analytical examples from published numerical experiments.
2. Map physical plant, feedback, predictor, target and assumptions to actual functions. A target-only integration is not a physical-plant reproduction.
3. Separate source parameters from repository choices, hypotheses, corrections and approximations. Record histories, horizons, steps, solvers, interpolation, guards and source revision.
4. Link measured CSV/JSON, source/environment records and exact test/CI evidence. Include all planned cases and failures. A test threshold is not automatically paper accuracy.
5. Explain finite-horizon outcomes and what is not established: asymptotic stability, theorem validity, certified region, robustness or exact figure agreement.
6. Separate original theory attribution from reproduction software. Preserve historical records in Git and explicitly supersede obsolete numerical claims.

Use a paper prefix for equation numbers when several sources are discussed. Label derived formulas; cite unnumbered paragraphs by page/section rather than inventing numbers. Refinement, replay of recorded commands, independent closed-loop checks and cross-language checks are different experiments. A red test remains a failure even when artifacts are saved; green CI applies only to executed checks at its tested source snapshot.
