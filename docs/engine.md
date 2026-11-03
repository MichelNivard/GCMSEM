# A general RAM-style cumulant engine

The derivative implementation is model-general. It does not contain separate formulas for a two-trait causal model, a factor model, or any other particular syntax. The parser translates a user's specification into path entries and sparse disturbance cumulant entries, each either fixed or mapped to a free parameter index. Several entries may share an index, implementing equality constraints.

Let `z` contain all observed and latent nodes and let `u` be their centered disturbances. The model is

\[
z=Az+u,\qquad y=Fz,\qquad T=F(I-A)^{-1}.
\]

`A` has outcome rows and predictor columns. Factor loadings are directed paths from a latent node to an indicator. `F` selects observed nodes. If `S` is the disturbance covariance matrix, the implied covariance is the usual RAM expression

\[
\Sigma=TST^\mathsf T.
\]

Let `Kq` be the symmetric disturbance cumulant tensor of order q. The same structure generalizes covariance to orders three and four:

\[
\kappa^{(q)}_{i_1\ldots i_q}
=\sum_{j_1,\ldots,j_q}
K^{(q)}_{j_1\ldots j_q}
\prod_{r=1}^{q}T_{i_rj_r}.
\]

Undeclared cross-disturbance cumulants are zero. An independent-source model therefore needs only marginal source cumulants. The sparse implementation also supports declared cross-source cumulants, with the distinct permutations of a symmetric tuple counted explicitly. For example, xxy has three distinct permutations, not one or six.

For every parameter, the derivative is

\[
\partial_\theta T
=F(I-A)^{-1}(\partial_\theta A)(I-A)^{-1},
\]

followed by

\[
\partial_\theta\kappa^{(q)}_I
=\sum_J(\partial_\theta K^{(q)}_J)\prod_rT_{i_rj_r}
+\sum_JK^{(q)}_J\sum_s
(\partial_\theta T_{i_sj_s})\prod_{r\ne s}T_{i_rj_r}.
\]

The product rule uses products excluding the differentiated factor; it never divides by a loading. Consequently exact zero loadings do not cause undefined derivatives. If a parameter labels multiple paths or source moments, all its contributions are accumulated. Setting a source's third and fourth cumulants to fixed Gaussian values removes those free parameters but uses the same engine.

This is a specialized exact first-derivative engine, rather than a general-purpose autodiff library. It supports user-defined **linear** SEM structures, including reciprocal paths if `I-A` is nonsingular. Stability and moment feasibility are additional diagnostics; matrix inversion alone does not establish either. Nonlinear structural equations, ordinal thresholds, arbitrary nonlinear parameter constraints, and mean structures are outside version 0.1's scope.

The native C++ backend accelerates the tensor contractions and their derivatives without an external tensor runtime. The R backend implements the same equations and serves as a reference for tests. Automatic parameter scaling improves optimization for very small SNP-effect moments while preserving model units and equality restrictions. The model Jacobian is reused for DWLS optimization, local identification checks and sandwich uncertainty.

For moment estimates `m`, diagonal inverse-variance weights `W`, and their complete sampling covariance `V`, the fit minimizes

\[
Q=(m-\mu(\theta))^\mathsf T W(m-\mu(\theta)),
\]

and uses

\[
\widehat{\mathrm{Var}}(\hat\theta)
=(J^\mathsf T WJ)^{-1}J^\mathsf T WVWJ(J^\mathsf T WJ)^{-1}.
\]

`V` is already the sampling covariance of the moment estimates; there is no further division by GWAS sample size.

The default `feasibility="check"` leaves DWLS unconstrained. Free negative disturbance variances and inadmissible higher cumulants are retained as classical or higher-order Heywood cases, with prominent warnings. The minimum-DWLS multistart solution is selected without using source admissibility as a tie-breaker or eligibility filter. Formal sandwich SEs remain available when numerically computable; source admissibility and inferential qualifications are separate flags. This distinguishes an inadmissible sample estimate from a failed optimizer or unidentified model, and avoids silently truncating the sampling distribution.

With explicit `feasibility="penalty"`, optimization adds an increasing quadratic penalty for violations of necessary source-moment PSD conditions. The penalty uses exact derivatives of the squared distance to the PSD cone, including at repeated eigenvalues, and fixed normalization for numerical scaling. Raw DWLS discrepancy and penalized objective are reported separately. This is an approximate feasibility strategy with finite penalty weights, not a proof of an exact global constrained optimum. Constraints change the estimator and can introduce boundary bias; they are not a default cleanup step. For large source systems, only covariance and marginal moment conditions are checked; the complete joint check is capped and its status is reported. Reciprocal stability remains an explicitly requested post-fit selection/check rather than a constraint parameterization.
