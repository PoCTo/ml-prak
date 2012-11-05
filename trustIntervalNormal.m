function X = trustIntervalNormal(P,sigma)
    X = norminv((P+1)/2,0,1)*sigma/sqrt(2);