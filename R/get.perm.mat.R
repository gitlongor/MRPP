nparts=function(n)
  as.bigz(factorialZ(sum(n))/prod(c(factorialZ(n), factorialZ(table(n)))))  ## total number of distinct trt assignments 

permuteTrt <-
function(trt, B=100L, idxOnly = FALSE) ## matrices of permutation vectors for one way design
{
    n=table(trt)
    cn=cumsum(n)
    ntrts=length(n)
    N=cn[ntrts]
    #mc=mchooseZ(N, n)
    ordn=order(n, names(n), decreasing=TRUE)
	trt=ordered(trt, levels=names(n)[ordn])

    SP=nparts(n)  
    part0=split(seq_len(N),trt); 
    
    if(B>=SP){ # list all partitions
        sp=setparts(n) # values are treatment indices
        B=ncol(sp)
#        for(i in seq(ntrts))  ans[[i]]=matrix( apply(sp==i,2L,function(xx)sort(which(xx)) ), ncol=B)
        ## the following 3 lines replace the previous line
			class(sp)=c('ordered','factor'); levels(sp)=names(n)[ordn]
			ans=split(rep.int(1:N,B),sp);  ## the previous line and this was originally implemented as ans=split(row(sp),sp) 
			for(i in seq(ntrts)) dim(ans[[i]])=c(length(ans[[i]])%/%B,B)   ## "/" returns double but %/% returns integer
        names(ans)=names(n)[ordn]
        
        #### swapping the original assignment to the first permutation        
        #flag=TRUE
        #for(b in seq(B)) if(setequal(part0, lapply(ans, '[', , b))){ idx=b; flag=FALSE; break }
        #if(isTRUE(flag)){
        #  warning("The first permutaiton may not be the original assignment.")
        #}
		#### the above block has be replaced by the following block for better speed
		part0.vec=do.call('c', part0); ans.vec=do.call('rbind',ans)
		b = which(.colSums(part0.vec == ans.vec, N, B) == N)
		stopifnot( length(b) ==1L )

        for(i in seq(ntrts)) {tmp=ans[[i]][,b]; ans[[i]][,b]=ans[[i]][,1L]; ans[[i]][,1L]=tmp} 

        if(isTRUE(idxOnly)){
			warning("'idxOnly=TRUE' has not been implemented yet when B is no larger than nparts(table(trt)). Full permutation vectors are returned.")		
        }#else{
            attr(ans, 'idx') = NA_character_
    		class(ans)='permutedTrt'
        #}
    }else{   #sample from all permutations using factoradic number. Ideally, a sample from 1:SP should work, but how to do this without enumerating all SP possibilities using setparts?

		facN = factorialZ(N)
		if(facN < as.bigz(6e15)){
				decfr=sample.bigz(facN, B) 	# faster for smaller facN
		}else   decfr=HSEL.bigz(facN, B)	# 	faster for larger facN
        idx=which(decfr==0L)
        if(length(idx)>0L) decfr[idx]=decfr[1L]
        decfr[1L]=as.bigz(0L)
		
		decfrCC=drop(as.character(decfr))  ## for speed only
		if(isTRUE(idxOnly)) {	## save memory
			ans = lapply(part0, as.matrix)
			attr(ans, 'idx') = decfrCC
			class(ans) = 'permutedTrt'
			return(ans)
		}
		
       ans=lapply(sapply(part0,length), matrix, data=NA_integer_, ncol=B)
		buff = integer(N); buff[1L]
       for(b in seq(B)){
            # perm=dec2permvec(decfr[b],N)  ## This subsetting decfr[b] is the slowest part!
            perm=dec2permvec(decfrCC[b],N)  ## This change speeds up for about 8~9X.
            # for(i in seq(ntrts)) ans[[i]][,b]=sort.int(perm[part0[[i]]])
			 for(i in seq(ntrts)) ans[[i]][,b]=.Call(radixSort_prealloc, perm[part0[[i]]], buff)  ## radix sort with pre-allocated buffer space
       }
       names(ans)=names(part0)
		attr(ans, 'idx') = NA_character_
		class(ans)='permutedTrt'
    }
    ans
}

nperms.permutedTrt=function(permutedTrt)
{
	if(is.na(attr(permutedTrt, 'idx')[1L])) {
			ncol(permutedTrt[[1L]])
	}else   length(attr(permutedTrt, 'idx'))
}

ntrt.permutedTrt=function(permutedTrt)
{
	sapply(permutedTrt, NROW)
}

trt.permutedTrt=function(permutedTrt)
{ ## previous implementation was bugged; this is the corrected version. 
	ans=rep(NA_integer_, sum(sapply(permutedTrt,NROW)))#, levels=seq_along(permutedTrt))
	for(i in seq_along(permutedTrt)) ans[permutedTrt[[i]][,1L]] = i
	class(ans)='factor'
	levels(ans)=names(permutedTrt)
	ans
}
