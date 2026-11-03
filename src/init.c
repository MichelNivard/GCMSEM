#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_gcmsem_contract(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"C_gcmsem_contract", (DL_FUNC) &C_gcmsem_contract, 6},
    {NULL, NULL, 0}
};
void attribute_visible R_init_GCMSEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
