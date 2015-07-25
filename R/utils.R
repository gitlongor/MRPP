safeseq =
function (from = 1L, to = 1L, by = 1L, ...) 
{
    disc = by * (from - to)
    if (disc > 0) {
        vector(class(disc), 0L)
    }
    else seq(from = from, to = to, by = by, ...)
}


midp = function(x, eps=1e-8)
{
	mean(x<x[1L]-eps) + .5* mean(x>=x[1L]-eps & x<=x[1L]+eps)
}
