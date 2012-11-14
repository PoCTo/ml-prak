function x = generateNormalRandomValueVonNeumanTrue(mu, sigma2)
    x=generateNormalRandomValueVonNeuman(mu,sigma2);
    while (x==Inf)
        x=generateNormalRandomValueVonNeuman(mu,sigma2);
    end
end
    