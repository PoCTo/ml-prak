function plotEmpiricalDistribution(x, opt)
    n = length(x);
    y = [0:n-1]./n;
    plot(sort(x),y, opt, 'LineWidth', 2);