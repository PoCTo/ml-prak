function plotEmpiricalExpectationBySize(x)
    n = length(x);
    xx = cumsum(x);
    %t = [1:n];
    xx = xx./[1:n];
    plot(1:n, xx);