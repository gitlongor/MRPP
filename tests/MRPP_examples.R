(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')
if(!R_CHECK_TIMINGS_){

if(interactive()){
#######  INSTALLATION OF MRPP PACKAGE VERSION 1.1-20170715-1120
setRepositories() # CRAN and BioC software should be selected at the minimum
chooseCRANmirror()
chooseBioCmirror()
## install prerequisites. On non-Windows systems, this might require manual installation of non-R software packages, e.g., gmp-devel, openssl-devel, etc. 
install.packages(c("gmp", "partitions", "QuasiSeq", "varComp", "moments", "PDQutils","devtools"))
## install MRPP package itself
devtools::install_github('gitlongor/MRPP', ref='20170715-1120')
}

library(MRPP)
#######  PERFORMING MRPP TESTS
set.seed(2340) #init random number seed
urand.bigz(0,seed=1032940L) # init big integer random number seed

n=10L; ntrt=2L; N=n*ntrt # 2 treatments, 10 indvidual in each
R=50L # number of variables
trt=gl(ntrt, n) # treatment assignment
x=matrix(rnorm(N*R),N) # simulated data: rows->individual; cols->variable
nparts( table(trt))		## 92378 partitions =  choose(20,10)/2

mrpp.obj=mrpp(x, trt, B=5e2L)	## creat an MRPP object, requesting only 500 random permutations
mrpp.test(mrpp.obj, method='permutation') # permutation p-value using 500 permutations

mrpp.obj=mrpp(x, trt, B=Inf)		## use all 92378 partitions 
mrpp.test(mrpp.obj, method='permutation')  ## truly exact test (slow)
mrpp.test(mrpp.obj, method='pearson3gca')  ## quick approx. p-value w/o permutations, fairly close. This can be used in modified MRPP with variable selection, since this approx. p-value is on a common scale, irrespective how many dimensions are selected. Very occasionally, the approx p-value might be outside [0,1].
mrpp.test(mrpp.obj, method='pearson3')  ## approx. p-value using only first 3 moments. This guarantees the p-value to be within [0,1]. 
mrpp.obj=mrpp(x, trt, B=100)		## use 100 perms
mrpp.test(mrpp.obj, method='pearson3gca.biweight')  ## this avoids out of bound p-values, by doing a small number of perms. This is a kind of semiparametri approximation to the permutation distribution. 


mrpp.test(x~trt, weight.trt='n') # formula interface with group weight being n
dist.x=dist(x)
permutedTrt=permuteTrt(trt, B=Inf) # slow part
mrpp.test(dist.x, permutedTrt=permutedTrt, method='permutation') # dist interface, exact test


####### CHOOSING BANDWIDTHS
permutedTrt=permuteTrt(trt, B=5000L) # could be slow if B is large; reuse this slow part
## 10 bandwidth selection method
# promising ones: 
bw.smoothp(x, permutedTrt, method='sym1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='dropadd1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='match.pear3', verbose=interactive())
bw.smoothp(x, permutedTrt, method='match.pear3gca', verbose=interactive())
bw.smoothp(x, permutedTrt, method='kde.mse1', verbose=interactive())
# others
bw.smoothp(x, permutedTrt, method='drop1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='add1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='dropaddsym1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='keep1', verbose=interactive())
bw.smoothp(x, permutedTrt, method='ss.gradp', verbose=interactive())

####### COMPUTING VARIABLE IMPORTANCE
## precompute to save time: 
mrpp.stats=mrpp.test(x,permutedTrt=permutedTrt,method='permutation')$all.statistics

## importance measures: 
iota=grad.smoothp(x, permutedTrt=permutedTrt, mrpp.stats=mrpp.stats, bw='match.pear3gca', adjust='none') # the bw argument can be any of the above 10 methods, or a precomputed numeric value
p.value(iota, type='drop1') # approx drop-1 p-value using computed derivatives
p.value(iota, type='add1') # approx add-1 p-value using computed derivatives
p.value(iota, type='keep1') # approx keep-1 p-value using computed derivatives

tau=grad.smoothp(x, permutedTrt=permutedTrt, mrpp.stats=mrpp.stats, bw=Inf, adjust='weighted.mean') # when bw=Inf, 'adjust' cannot be 'none'

cor(iota, tau); cor(iota, tau, method='spear') # compare the two

iota_approx_to_tau=grad.smoothp(x, permutedTrt=permutedTrt, mrpp.stats=mrpp.stats, bw=100, adjust='weighted.mean') # large but finite bw. When adjust is weighted.mean, iota and tau are on similar scales. 
plot(iota_approx_to_tau, tau); abline(0,1) # compare


}
