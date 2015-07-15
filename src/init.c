#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern SEXP radixSort_prealloc(SEXP, SEXP); 
extern SEXP radixSort1PassByCol(SEXP, SEXP); 
extern SEXP mrppstats(SEXP , SEXP , SEXP ); 
extern SEXP FR2permvec(const SEXP , const SEXP);

R_CallMethodDef callMethods[]  = {
       {"radixSort1PassByCol", (DL_FUNC) &radixSort1PassByCol, 2},
       {"radixSort_prealloc", (DL_FUNC) &radixSort_prealloc, 2},
       {"mrppstats", (DL_FUNC) &mrppstats, 3},
       {"C_FR2permvec", (DL_FUNC) &FR2permvec, 2},
       {NULL, NULL, 0}
};

void R_init_MRPP(DllInfo *info)
 {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
 }
 