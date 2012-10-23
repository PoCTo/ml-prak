function x = generateKantorRandomSequence(prec, sz)
    x = zeros(1, sz);
    for i = 1:sz 
        x(i) = generateKantorRandomValue(prec);
    end