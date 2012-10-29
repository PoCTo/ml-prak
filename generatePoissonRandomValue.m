function x = generatePoissonRandomValue(l)
    sum = 1;
    i = 0;
    res = exp(-l)
    while (sum >= res)
        sum = sum * rand(1)
        i = i+1;
    end
    x = i-1;