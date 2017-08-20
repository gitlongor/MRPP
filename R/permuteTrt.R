nparts=function(n)
  as.bigz(factorialZ(sum(n))/prod(c(factorialZ(n), factorialZ(table(n)))))  ## total number of distinct trt assignments 

eqv.trt=function(trt)
{
	trt.bak=trt
	trt=unclass(trt)

	n=table(trt)
	n.names=as.integer(names(n))
	
	ntrt=length(n)
	N=as.integer(sum(n))
	tabn=table(n)
	tabn.names=as.integer(names(tabn))
	
	eqv.trts=matrix(trt, N, prod(factorial(tabn)))
	
	neqv.trts=NCOL(eqv.trts)
	nperms=1L
	for(i in seq_along(tabn)){
		if(tabn[i]==1L) next
		idx=which(n==tabn.names[i])
		perm.trt=idx[perms(tabn[i])]
		dim(perm.trt)=c(tabn[i], factorial(tabn[i]))
		for(j in seq_len(NCOL(perm.trt))){
			for(k in seq_len(tabn[i]))
				eqv.trts[trt==idx[k], nperms]=perm.trt[k,j]
			nperms=nperms+1L
		}
	}
	stopifnot(sum(duplicated.matrix(apply(eqv.trts,2L,table),MARGIN=2L))==neqv.trts-1L)
	levels(eqv.trts)=levels(trt.bak)
	class(eqv.trts)=class(trt.bak)
	eqv.trts
}

permuteTrt <-
function(trt, B=100L, idxOnly = FALSE, sample.method=c('sample','permute')) ## matrices of permutation vectors for one way design
{
	sample.method=match.arg(sample.method)
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
		## locating/swapping the original treatment assignment
		b=seq_len(B)
		for(i in seq_along(n)){ ## pretty fast as the size of b is keeping shrinking
			if(n[i]<=1L)break;
			j0=part0[[i]][1]
			for(j in part0[[i]][-1L])
				b=b[which(sp[j0,b]==sp[j,b])]
		}
		stopifnot(length(b)==1L)
		tmp=sp[,b]; sp[,b]=sp[,1L]; sp[,1L]=tmp
		
		## matching sp treatment indices with treatment labels
		levs=levels(trt)
		levs[sp[sapply(part0,'[[',1L),1L]]=levs
		class(sp)='factor'; attr(sp,'levels')=levs
		ans=split(rep.int(1:N,B),sp);  
		## the previous block was originally implemented as ans=split(row(sp),sp) 
		for(i in seq_len(ntrts)) dim(ans[[i]])=c(length(ans[[i]])%/%B,B)   ## "/" returns double but %/% returns integer
		## the previous block was initially implemented as for(i in seq(ntrts))  ans[[i]]=matrix( apply(sp==i,2L,function(xx)sort(which(xx)) ), ncol=B)
		
		ans = ans[levels(trt)]
		stopifnot(identical(part0, lapply(ans,'[',,1L)))
			
        if(isTRUE(idxOnly)){
			warning("'idxOnly=TRUE' has not been implemented yet when B is no larger than nparts(table(trt)). Full permutation vectors are returned.")		
        }#else{
            attr(ans, 'idx') = NA_character_
    		class(ans)='permutedTrt'
        #}
    }else if(sample.method=='sample')
	{   #sample from all permutations using factoradic number. Ideally, a sample from 1:SP should work, but how to do this without enumerating all SP possibilities using setparts?

		facN = factorialZ(N)
		if(facN < as.bigz(6e15)){
				decfr=sample.bigz(facN, B) 	# faster for smaller facN
		}else   decfr=HSEL.bigz(facN, B)	# 	faster for larger facN
		idx1=which(decfr==0L)
		
		decfrCC = as.character(decfr)
		if(isTRUE(idxOnly)) {	## save memory
			ans = lapply(part0, as.matrix)
				if(length(idx1)>0L) decfrCC[idx1]=decfrCC[1L]
				decfrCC[1L]='0'
			attr(ans, 'idx') = decfrCC
			class(ans) = 'permutedTrt'
			return(ans)
		}
		
		perms = sapply(decfrCC, dec2permvec, N=N)
			if(length(idx1)>0L) perms[,idx1]=perms[,1L]
			perms[,1L] = seq_len(N)
    }else if(sample.method=='permute'){
		if(FALSE){ #implementation based on the hashing of environemnts
			perm.env=new.env(hash=TRUE, parent=globalenv(),size=B)
			nperm.found=0L
			char.lab=as.character(seq_len(N))
			repeat{
				tmp=paste(sample(char.lab),collapse=',')
				if(!is.null(perm.env[[tmp]])) next
				delayedAssign(tmp, TRUE, assign.env=perm.env)
				nperm.found=nperm.found+1L
				if(nperm.found==B) break
			}
			
			# alternative implementation
			char.lab=as.character(seq_len(N))
			perm.env=new.env(hash=TRUE, parent=globalenv(),size=B)
			n.idx=B-1L
			delayedAssign(paste(seq_len(N),collapse=','),TRUE, assign.env=perm.env)
			repeat{
				tmp=replicate(n.idx, paste(sample(char.lab),collapse=','))
				for(i in seq_len(n.idx))delayedAssign(tmp[i],TRUE, assign.env=perm.env)
				n.idx=B-length(ls(perm.env))
				if(n.idx==0) break
			}
			
		}
		if(TRUE){ #implementation based on hashing of uniqueAtomMat
			tabn=table(n)
			eqv.trts=eqv.trt(trt)
			neqv.trts=NCOL(eqv.trts)
			
			tmp.idx=2:B
			n.idx=B-1L
			dim.bak=c(N,B)
			ans=array(NA_integer_, dim=dim.bak)
			ans[,1L]=seq_len(N)
			repeat{
				ans[,tmp.idx]=replicate(n.idx, sample.int(N))
				dim(ans)=NULL
				tmp=eqv.trts[ans,]
				dim(tmp)=c(N, neqv.trts*B)
				tmp.idx=which(uniqueAtomMat::duplicated.matrix(tmp,MARGIN=2L))
				dim(ans)=dim.bak
				tmp.idx=unique(tmp.idx%%B+1L)
				tmp.idx=tmp.idx[tmp.idx!=1L]
				if(length(tmp.idx)==0) break
				n.idx=length(tmp.idx)
			}
		}
		
	}else stop('unknown "sample.method"')

	ans = split.data.frame(perms, trt)
   for(i in seq_len(ntrts)) ans[[i]] = .Call(radixSort1PassByCol, ans[[i]], N)
   
	attr(ans, 'idx') = NA_character_
	class(ans)='permutedTrt'
	
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
{
	ans=rep(NA_integer_, sum(sapply(permutedTrt,NROW)))#, levels=seq_along(permutedTrt))
	for(i in seq_along(permutedTrt)) ans[permutedTrt[[i]][,1L]] = i
	class(ans)='factor'
	attr(ans, 'levels')=names(permutedTrt)
	ans
}
