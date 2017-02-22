#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP causalInference(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dagToEssentialGraph(SEXP, SEXP);
extern SEXP estimateSkeleton(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP globalMLE(SEXP, SEXP, SEXP, SEXP);
extern SEXP globalScore(SEXP, SEXP, SEXP, SEXP);
extern SEXP localMLE(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP localScore(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optimalTarget(SEXP, SEXP);
extern SEXP representative(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"causalInference",     (DL_FUNC) &causalInference,     5},
    {"dagToEssentialGraph", (DL_FUNC) &dagToEssentialGraph, 2},
    {"estimateSkeleton",    (DL_FUNC) &estimateSkeleton,    7},
    {"globalMLE",           (DL_FUNC) &globalMLE,           4},
    {"globalScore",         (DL_FUNC) &globalScore,         4},
    {"localMLE",            (DL_FUNC) &localMLE,            5},
    {"localScore",          (DL_FUNC) &localScore,          5},
    {"optimalTarget",       (DL_FUNC) &optimalTarget,       2},
    {"representative",      (DL_FUNC) &representative,      1},
    {NULL, NULL, 0}
};

void R_init_pcalg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
