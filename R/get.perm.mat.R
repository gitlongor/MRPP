nparts=function(n)
  as.bigz(factorialZ(sum(n))/prod(c(factorialZ(n), factorialZ(table(n)))))  ## total number of distinct trt assignments 

permuteTrt <-
function(trt, B=100L, idxOnly = FALSE) ## matrices of permutation vectors for one way design
{
    n=table(trt)
    cn=cumsum(n)
    ntrts=length(n)
    N=cn[ntrts]
    ordn=order(-n, names(n), decreasing=FALSE)
	trt=ordered(trt, levels=names(n)[ordn])
	trtc=as.character(trt)
	n=n[ordn]

    SP=nparts(n)  
    part0=split(seq_len(N),trt); 
    
    if(B>=SP){ # list all partitions
        sp=setparts(n) # values are treatment indices
        B=ncol(sp)
			## matching sp group labels with desired label orders
            b=seq(B)
            for(i in seq_along(n)){
                if(n[i]<=1L)break;
                j0=part0[[i]][1]
                for(j in part0[[i]][-1L])
                    b=b[which(sp[j0,b]==sp[j,b])]
            }
            stopifnot(length(b)==1L)
            tmp=sp[,b]; sp[,b]=sp[,1L]; sp[,1L]=tmp
            
#        for(i in seq(ntrts))  ans[[i]]=matrix( apply(sp==i,2L,function(xx)sort(which(xx)) ), ncol=B)
        ## the following 3 lines replace the previous line
            levs=levels(trt)
            levs[sp[sapply(part0,'[[',1L),1]]=levs
            class(sp)='factor'; attr(sp,'levels')=levs
			ans=split(rep.int(1:N,B),sp);  ## the previous line and this was originally implemented as ans=split(row(sp),sp) 
			for(i in seq(ntrts)) dim(ans[[i]])=c(length(ans[[i]])%/%B,B)   ## "/" returns double but %/% returns integer
			ans = ans[levels(trt)]
			stopifnot(identical(part0, lapply(ans,'[',,1L)))
			
		if(FALSE){
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
		}

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
		idx1=which(decfr==0L)
		
		if(isTRUE(idxOnly)) {	## save memory
			ans = lapply(part0, as.matrix)
			decfrCC = as.character(decfr)
				if(length(idx1)>0L) decfrCC[idx1]=decfrCC[1L]
				decfrCC[1L]='0'
			attr(ans, 'idx') = decfrCC
			class(ans) = 'permutedTrt'
			return(ans)
		}
		
		perms = sapply(decfr, dec2permvec, N=N)
			if(length(idx1)>0L) perms[,idx1]=perms[,1L]
			perms[,1L] = seq(N)
		ans = split.data.frame(perms, trt)
	   for(i in seq(ntrts)) ans[[i]] = .Call(radixSort1PassByCol, ans[[i]], N)
	   
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
