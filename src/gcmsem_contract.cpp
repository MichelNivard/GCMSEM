// Sparse symmetric cumulant contraction and exact first derivatives.
// No finite differences and no divisions by loadings: zero paths remain regular.
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <vector>

namespace {
struct Source {
    double value;
    int parameter; // zero based; -1 denotes fixed
    std::vector<std::array<int, 4> > permutations;
};
SEXP member(SEXP object, const char* name) {
    SEXP names = Rf_getAttrib(object, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(object); ++i)
        if (names != R_NilValue && !strcmp(CHAR(STRING_ELT(names, i)), name))
            return VECTOR_ELT(object, i);
    throw std::invalid_argument("Missing source-row field");
}
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
SEXP contract(SEXP A_, SEXP dA_, SEXP K_, SEXP tuples_, SEXP p_, SEXP jacobian_) {
    require(TYPEOF(A_) == REALSXP && Rf_isMatrix(A_), "A must be a double matrix");
    SEXP dims = Rf_getAttrib(A_, R_DimSymbol);
    const int nobs = INTEGER(dims)[0], nsource = INTEGER(dims)[1];
    require(nobs > 0 && nsource > 0, "A must have positive dimensions");
    const int p = Rf_asInteger(p_);
    const int deriv = Rf_asLogical(jacobian_);
    require(p != NA_INTEGER && p >= 0 && deriv != NA_LOGICAL, "Invalid parameter count or jacobian flag");
    require(TYPEOF(K_) == VECSXP && XLENGTH(K_) == 4, "K must have four order-indexed entries");
    require(TYPEOF(tuples_) == VECSXP, "tuples must be a list");
    const R_xlen_t nmoment = XLENGTH(tuples_);
    require(nmoment <= INT_MAX, "Too many requested moments for an R matrix");
    const double* A = REAL(A_);
    for (R_xlen_t i = 0; i < XLENGTH(A_); ++i)
        require(R_FINITE(A[i]), "Non-finite entry in A");
    const double* dA = NULL;
    const R_xlen_t plane = static_cast<R_xlen_t>(nobs) * nsource;
    if (deriv) {
        SEXP dd = Rf_getAttrib(dA_, R_DimSymbol);
        require(TYPEOF(dA_) == REALSXP && XLENGTH(dd) == 3 &&
                INTEGER(dd)[0] == nobs && INTEGER(dd)[1] == nsource && INTEGER(dd)[2] == p,
                "dA dimensions must be nrow(A), ncol(A), p");
        dA = REAL(dA_);
        for (R_xlen_t i = 0; i < XLENGTH(dA_); ++i)
            require(R_FINITE(dA[i]), "Non-finite entry in dA");
    }
    std::array<std::vector<Source>, 5> sources;
    for (int q = 2; q <= 4; ++q) {
        SEXP rows = VECTOR_ELT(K_, q - 1);
        require(rows == R_NilValue || TYPEOF(rows) == VECSXP, "K order entries must be lists");
        for (R_xlen_t r = 0; r < XLENGTH(rows); ++r) {
            SEXP row = VECTOR_ELT(rows, r);
            require(TYPEOF(row) == VECSXP, "Source rows must be lists");
            Source source;
            source.value = Rf_asReal(member(row, "value"));
            const int parameter = Rf_asInteger(member(row, "parameter"));
            require(R_FINITE(source.value), "Source cumulants must be finite");
            require(parameter != NA_INTEGER && parameter >= 0 && parameter <= p, "Source parameter index out of range");
            source.parameter = parameter - 1;
            if (source.value == 0 && (!deriv || parameter == 0)) continue;
            SEXP tuple = member(row, "tuple");
            require(TYPEOF(tuple) == INTSXP && XLENGTH(tuple) == q, "Source tuple must have the correct integer order");
            std::array<int, 4> perm = {{0, 0, 0, 0}};
            for (int k = 0; k < q; ++k) {
                const int v = INTEGER(tuple)[k];
                require(v != NA_INTEGER && v >= 1 && v <= nsource, "Source tuple index out of range");
                perm[k] = v - 1;
            }
            std::sort(perm.begin(), perm.begin() + q);
            do { source.permutations.push_back(perm); }
            while (std::next_permutation(perm.begin(), perm.begin() + q));
            sources[q].push_back(source);
        }
    }
    std::vector<std::array<int, 4> > observed(nmoment);
    std::vector<int> orders(nmoment);
    for (R_xlen_t m = 0; m < nmoment; ++m) {
        SEXP tuple = VECTOR_ELT(tuples_, m);
        const int q = static_cast<int>(XLENGTH(tuple));
        require(TYPEOF(tuple) == INTSXP && q >= 2 && q <= 4, "Observed tuples must have integer order two through four");
        orders[m] = q;
        for (int k = 0; k < q; ++k) {
            const int v = INTEGER(tuple)[k];
            require(v != NA_INTEGER && v >= 1 && v <= nobs, "Observed tuple index out of range");
            observed[m][k] = v - 1;
        }
    }
    SEXP values = PROTECT(Rf_allocVector(REALSXP, nmoment));
    SEXP jac = PROTECT(deriv ? Rf_allocMatrix(REALSXP, static_cast<int>(nmoment), p) : R_NilValue);
    std::fill(REAL(values), REAL(values) + nmoment, 0.0);
    double* J = deriv ? REAL(jac) : NULL;
    if (deriv) std::fill(J, J + nmoment * p, 0.0);
    for (R_xlen_t m = 0; m < nmoment; ++m) {
        if (m % 256 == 0) R_CheckUserInterrupt();
        const int q = orders[m];
        for (const Source& source : sources[q]) {
            double coefficient = 0;
            for (const auto& perm : source.permutations) {
                std::array<double, 4> a;
                std::array<R_xlen_t, 4> offsets;
                double product = 1;
                for (int k = 0; k < q; ++k) {
                    offsets[k] = observed[m][k] + static_cast<R_xlen_t>(nobs) * perm[k];
                    a[k] = A[offsets[k]];
                    product *= a[k];
                }
                coefficient += product;
                if (deriv && source.value != 0) {
                    for (int k = 0; k < q; ++k) {
                        double weight = source.value;
                        for (int other = 0; other < q; ++other)
                            if (other != k) weight *= a[other];
                        if (weight == 0) continue;
                        for (int j = 0; j < p; ++j)
                            J[m + nmoment * j] += weight * dA[offsets[k] + plane * j];
                    }
                }
            }
            REAL(values)[m] += source.value * coefficient;
            if (deriv && source.parameter >= 0)
                J[m + nmoment * source.parameter] += coefficient;
        }
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, values);
    SET_VECTOR_ELT(result, 1, jac);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("estimate"));
    SET_STRING_ELT(names, 1, Rf_mkChar("jacobian"));
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(4);
    return result;
}
} // namespace

extern "C" SEXP C_gcmsem_contract(SEXP A, SEXP dA, SEXP K, SEXP tuples, SEXP p, SEXP jacobian) {
    try {
        return contract(A, dA, K, tuples, p, jacobian);
    } catch (const std::exception& e) {
        Rf_error("GCMSEM native contraction: %s", e.what());
    } catch (...) {
        Rf_error("GCMSEM native contraction: unknown C++ exception");
    }
    return R_NilValue;
}
