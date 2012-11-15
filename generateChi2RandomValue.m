function x = generateChi2RandomValue(k)
    r = rand(1,2);
    r1=r(1); r2=r(2);
    seq1 = arrayfun(@(x)generateNormalRandomValuePair1(0,1,r1,r2),1:k);
    seq1 = seq1.*seq1;
    x = sum(seq1);
    