#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern SEXP radixSort_prealloc(SEXP, SEXP); 
extern SEXP radixSort1PassByCol(SEXP, SEXP); 
extern SEXP radixSort_preallocMax(SEXP, SEXP, SEXP); 
extern SEXP mrppstats(SEXP , SEXP , SEXP ); 
extern SEXP FR2permvec(const SEXP , const SEXP);
extern SEXP pmax0(SEXP);


R_CallMethodDef callMethods[]  = {
       {"radixSort1PassByCol", (DL_FUNC) &radixSort1PassByCol, 2},
       {"radixSort_prealloc", (DL_FUNC) &radixSort_prealloc, 2},
       {"radixSort_preallocMax", (DL_FUNC) &radixSort_preallocMax, 3},
       {"mrppstats", (DL_FUNC) &mrppstats, 3},
       {"C_FR2permvec", (DL_FUNC) &FR2permvec, 2},
       {"pmax0", (DL_FUNC) &pmax0, 1},
       {NULL, NULL, 0}
};

void R_init_MRPP(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
