function X=generateParetoRandomMatrix(x_m,k,sz_x,sz_y)
    X = rand(sz_x, sz_y);
    X = X.^(-1/k)*x_m;
    
    