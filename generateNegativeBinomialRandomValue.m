function x = generateNegativeBinomialRandomValue(r, p)
    count = 0;
    index = 0;
    while (count < r)
        x = generateBernulliRandomValue(p);
        index = index + 1;
        if (x == 1)
            count = count + 1;
        end
    end
    x = index - count;