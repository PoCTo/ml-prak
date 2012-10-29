function bars = generateRightNegativeBinomialCount(r, p, sz)
bars=[]
i=0;
res = 1;
eps = 0.001;

pr=p^r;
qk=1;
while res > eps
    res = nchoosek(i+r-1, i)*qk*pr;
    qk=qk*(1-p);
    bars=[bars, res];
    i=i+1;
end
bars = bars.*sz;
