function x = generateNormalRandomValuePair(mu, sigma2)
    y = rand(1,2);
    x = sigma2*sqrt(-2*log(y(1)))*...
        [sin(2*pi*y(2));cos(2*pi*y(2))]...
        +mu;