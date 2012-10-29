function x = generateBinomialRandomSequence(n, p, sz)
    x = zeros(1, sz);
    for i = 1:sz
        x(i) = generateBinomialRandomValue(n, p);
    end