function x = generateChi2RandomValuePair1(k,r1,r2)
    seq1 = arrayfun(@(x)generateNormalRandomValuePair1(0,1,r1,r2),1:k);
    seq1 = seq1.*seq1;
    x = sum(seq1);
    