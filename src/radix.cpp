#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define RADIXEXP (3) 
template <int n>
struct radix {
	enum { value = 2 * radix<n - 1>::value };
};
template <>
struct radix<1> {
	enum { value = 2 };
};


extern "C" void radixsort(int *a, int n, int *b)
{
	static int bucket[radix<RADIXEXP>::value]; 
	const size_t bucketSize = sizeof(int) * radix<RADIXEXP>::value;
	int m = a[n-1];
	int * tmp;
	for (tmp = a + n -2; tmp>=a; --tmp) if (*tmp > m) m = *tmp;

	const unsigned int FFs = radix<RADIXEXP>::value -1;
	bool flag = false;
	int * a0 = a;
	int * abends[2]={a+n, b+n};
	int * aend = abends[0];
	int * buckEnd = bucket + radix<RADIXEXP>::value;
	int * current;
	for(int shift=0 ; (m >> shift) > 0 ; shift += RADIXEXP)
	{
		memset(bucket, 0, bucketSize);
		for (tmp = aend-1; tmp >= a; --tmp)
			++bucket[(*tmp >> shift) & FFs];
		for (tmp=bucket,current=tmp+1; current != buckEnd; ++tmp,++current)
			*current += *tmp;
		for (tmp = aend-1; tmp >= a; --tmp)
			b[--bucket[(*tmp >> shift) & FFs]] = *tmp;
		//for (i = 0; i < n; i++)      a[i] = b[i]; 
		tmp=a; a=b; b=tmp; //swap a and b pointers to avoid copy data
		flag = !flag;
		aend = abends[flag];
	}
	
	if(a!=a0) memcpy(a0, a, sizeof(int)*n);
}

void radixsort(int *a, int n, int *b, int m) // with known(or approx) maximum
{
	static int bucket[radix<RADIXEXP>::value]; 
	const size_t bucketSize = sizeof(int) * radix<RADIXEXP>::value;
	int * tmp;
	//for (tmp = a + n -2; tmp>=a; --tmp) if (*tmp > m) m = *tmp;

	const unsigned int FFs = radix<RADIXEXP>::value -1;
	bool flag = false;
	int * a0 = a;
	int * abends[2]={a+n, b+n};
	int * aend = abends[0];
	int * buckEnd = bucket + radix<RADIXEXP>::value;
	int * current;
	for(int shift=0 ; (m >> shift) > 0 ; shift += RADIXEXP)
	{
		memset(bucket, 0, bucketSize);
		for (tmp = aend-1; tmp >= a; --tmp)
			++bucket[(*tmp >> shift) & FFs];
		for (tmp=bucket,current=tmp+1; current != buckEnd; ++tmp,++current)
			*current += *tmp;
		for (tmp = aend-1; tmp >= a; --tmp)
			b[--bucket[(*tmp >> shift) & FFs]] = *tmp;
		//for (i = 0; i < n; i++)      a[i] = b[i]; 
		tmp=a; a=b; b=tmp; //swap a and b pointers to avoid copy data
		flag = !flag;
		aend = abends[flag];
	}
	
	if(a!=a0) memcpy(a0, a, sizeof(int)*n);
}

extern "C" SEXP radixSort_prealloc(SEXP x, SEXP buff)
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

extern "C" SEXP radixSort_preallocMax(SEXP x, SEXP buff, SEXP m)
// Radix sort using pre-allocated buffer space. 
// m is is assumed to be >= max(x)
// For speed, this does not check any assertions.
{
	int * ptrAns;
	R_len_t n;
	SEXP ans;
	
	n=LENGTH(x);
	
	PROTECT(ans = NEW_INTEGER(n));
	ptrAns = INTEGER(ans);
	memcpy (ptrAns, INTEGER(x), sizeof(int) * n);
	
	radixsort(ptrAns, n, INTEGER(buff), *INTEGER(m));
	
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

extern "C" SEXP radixSort1PassByCol(SEXP x, SEXP maxx)
/* Single-pass radix sort using a common buffer space for each col of x.
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
