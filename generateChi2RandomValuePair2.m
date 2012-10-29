function x = generateChi2RandomValuePair2(k)
    seq1 = arrayfun(@(x)generateNormalRandomValuePair2(0,1),1:k);
    seq1 = seq1.*seq1;
    x = sum(seq1);
    