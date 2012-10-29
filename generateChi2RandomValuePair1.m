function x = generateChi2RandomValuePair1(k)
    seq1 = arrayfun(@(x)generateNormalRandomValuePair1(0,1),1:k);
    seq1 = seq1.*seq1;
    x = sum(seq1);
    