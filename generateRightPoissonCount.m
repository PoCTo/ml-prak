function x = generateRightPoissonCount(l,sz)
epp = 1;
x = [];
i = 0;
fac = 1;
li=1;
while i<2*l+2
    epp = li/fac*exp(-l)
    x = [x, epp*sz];
    i=i+1;
    fac = fac*i;
    li=li*l;
end