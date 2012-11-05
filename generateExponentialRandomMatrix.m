function x = generateExponentialRandomMatrix(l,sz_x,sz_y)
    x = rand(sz_x,sz_y);
    x = - log(x)/l;