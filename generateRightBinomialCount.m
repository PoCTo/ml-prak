function bars = generateRightBinomialCount(n, p, sz)
bars = arrayfun(@(x)fnc(n,x,p),0:n).*sz;

function res = fnc(n,k,p)
res = nchoosek(n,k)*p^k*(1-p)^(n-k);
