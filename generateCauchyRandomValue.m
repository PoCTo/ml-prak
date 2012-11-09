    function x = generateCauchyRandomValue(x0, gamma)
    x = rand(1);
    x = x0+gamma*tan(pi*(x-1/2));