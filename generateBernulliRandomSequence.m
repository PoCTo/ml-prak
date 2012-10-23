% p - вероятность единицы
% n - количество величин
function x = generateBernulliRandomSequence(p, sz)
    x = rand(1, sz);
    t = find(x > 1-p);
    x(t) = ones(size(t));
    t = find(x <= 1-p);    
    x(t) = zeros(size(t));