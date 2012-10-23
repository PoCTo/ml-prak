function x = generateBernulliRandomValue(p)
    x = rand(1,1);
    if (x > p)
        x = 0;
    else
        x = 1;
    end