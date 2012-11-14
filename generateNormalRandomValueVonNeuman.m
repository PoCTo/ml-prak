function x = generateNormalRandomValueVonNeuman(mu, sigma2)
    k = sqrt(2*pi/exp(1));
    x = generateCauchyRandomValue(0,1);
    check = generateBernulliRandomValue(sqrt(pi/2)*exp(-x^2/2)*...
        (1+x^2)/k);
    if (check == 0) 
        x=Inf;
    else
        x = mu+sigma2*x;
    end