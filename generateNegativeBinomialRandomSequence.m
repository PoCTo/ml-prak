function x = generateNegativeBinomialRandomSequence(r, p, sz)
    x = zeros(1, sz);
    for i = 1:sz 
        x(i) = generateNegativeBinomialRandomValue(r, p);
    end