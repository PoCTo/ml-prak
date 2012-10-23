% prec - precision
function x = generateKantorRandomValue(prec)
    x = 0;
    multiplier = 1/3;
    for i = 1:prec
        t = generateBernulliRandomValue(0.5);
        if (t == 1)
            t = 2;
        end
        x = x + t * multiplier;
        multiplier = multiplier / 3;
    end