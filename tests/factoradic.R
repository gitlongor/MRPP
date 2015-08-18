(R_CHECK_TIMINGS_ = Sys.getenv('_R_CHECK_TIMINGS_') != '')


library(MRPP)
stopifnot(FR2dec(dec2FR(10,10))==10)
stopifnot(permvec2dec(dec2permvec(100,10))==100)

set.seed(2340)
pv=sample(10L)
stopifnot(all(FR2permvec(permvec2FR(pv))==pv))
