# Contributing

Keep the numerical estimator, reference-score construction and SEM fit separately testable. Changes to the estimand, effect scale, SNP orientation, moment ordering, or jackknife covariance require a schema or version change and a migration note.

Run Python unit tests and R package checks before proposing a change. Add tests that demonstrate a scientific identity, independently calculated expected value, malformed input rejection, or a previously observed defect. Compare exact R Jacobians with finite differences in tests, not during fitting.

Numerical changes should retain deterministic small fixtures and run the GWASBrewer benchmark. Record failures, boundaries, and coverage rather than selecting favorable seeds. A faster optimizer is useful only if it preserves model predictions and derivative correctness.

Use issues to report the software version, command or model, relevant log lines, operating system, and a minimal synthetic example. Do not upload participant data, restricted summary statistics, credentials, or licensed reference panels.

