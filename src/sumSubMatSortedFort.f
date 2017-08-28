      double precision function sumSubMatSortedFort(x, idx, n, nsamp)
      integer idx(*), n, nsamp
      double precision x(*) 

c     Local variables
	  integer i, j, n2, base
	  
	  sumSubMatSortedFort = 0.0
	  n2 = 2*nsamp 

	  do 200 j=1, n-1
c       column index shifted by idx(j)
	    base =  (n2-idx(j)) * (idx(j) -1) / 2 - idx(j)
		do 100 i=j+1, n
c          Direct summation algorithm
          sumSubMatSortedFort = sumSubMatSortedFort +x(base+idx(i))
 100    continue		
 200  continue
      sumSubMatSortedFort = sumSubMatSortedFort * 2.0
      return
      end


      double precision function sumSubMatSrtFtKahan(x, idx, n, nsamp)
      integer idx(*), n, nsamp
      double precision x(*) 

c     Local variables
      double precision compen, adjx, tmpAns
	  integer i, j, n2, base
	  
	  sumSubMatSrtFtKahan = 0.0
	  compen = 0.0
	  n2 = 2*nsamp 

	  
	  do 1200 j=1, n-1
c       column index shifted by idx(j)
	    base =  (n2-idx(j)) * (idx(j) -1) / 2 - idx(j)
		do 1100 i=j+1, n
c         Kahan's summation algorithm
		  adjx = x(base + idx(i)) - compen
		  tmpAns = sumSubMatSrtFtKahan +adjx 
		  compen = tmpAns - sumSubMatSrtFtKahan - adjx
		  sumSubMatSrtFtKahan = tmpAns
 1100   continue		
 1200 continue
      sumSubMatSrtFtKahan = sumSubMatSrtFtKahan * 2.0
      return
      end
