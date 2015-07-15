#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define RADIX (8) 

void radixsort(int *a, int n, int *b)
{
  static int bucket[RADIX]; 
  int i, m = a[0], exp = 1;
  for (i = 0; i < n; ++i)
  {
    if (a[i] > m)
      m = a[i];
  }
  
  
  while (m / exp > 0)
  {
	for(i=0; i<RADIX; ++i) bucket[i]=0;  // memset (bucket, 0, RADIX); 
    // int bucket[RADIX] =  {  0 };
    for (i = 0; i < n; ++i)
      bucket[a[i] / exp % RADIX]++;
    for (i = 1; i < RADIX; ++i)
      bucket[i] += bucket[i - 1];
    for (i = n - 1; i >= 0; --i)
      b[--bucket[a[i] / exp % RADIX]] = a[i];
    for (i = 0; i < n; i++)
      a[i] = b[i];
    exp *= RADIX;
   }
}

SEXP radixSort_prealloc(SEXP x, SEXP buff)
// Radix sort using pre-allocated buffer space. 
// For speed, this does not check any assertions.
{
	int * ptrAns;
	R_len_t n;
	SEXP ans;
	
	n=LENGTH(x);
	
	PROTECT(ans = NEW_INTEGER(n));
	ptrAns = INTEGER(ans);
	memcpy (ptrAns, INTEGER(x), sizeof(int) * n);
	
	radixsort(ptrAns, n, INTEGER(buff));
	
	UNPROTECT(1);
	return ans;
}


static inline void radixsort1pass(int *a, int n, int max1p, int *out, int* buff)
{
	int i;
  
	memset (buff, 0, max1p * sizeof(int) ); // buff is one-based
	for (i = n-1; i >= 0; --i)	  buff[a[i]]++;
	for (i = 2; i < max1p; ++i)	  buff[i] += buff[i - 1];
	for (i = n-1; i >= 0; --i)	  out[--buff[a[i]]] = a[i];
}

SEXP radixSort1PassByCol(SEXP x, SEXP maxx)
/* Single-pass radix sort using pre-allocated buffer space for each col of x.
 * For speed, this does not check any assertions.
 * This is only used when inputs are relatively small integers in the range of 1..maxx
 */ 
{
	int * pIn, * pOut, *pBuff, * out0;
	R_len_t n;
	SEXP out, buff;
	int* dim;
	int m, shift;
	
	m = INTEGER(maxx)[0]+1; // buff is one-based
	dim=INTEGER(getAttrib(x, R_DimSymbol));	
	n = dim[0];
	out = PROTECT(allocMatrix(INTSXP, dim[0], dim[1]));
	buff = PROTECT(allocVector(INTSXP, m));
	pBuff = INTEGER(buff);
	
	shift = (dim[1] - 1) * n;
	pIn = INTEGER(x) + shift;
	out0 = INTEGER(out);
	for (pOut = out0 + shift; pOut >= out0; pOut -= n){
		radixsort1pass(pIn, n, m, pOut, pBuff);
		pIn -= n; 
	}
	
	UNPROTECT(2);
	return out;
}
