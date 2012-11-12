%% 7 (eps,p) -> n
p = 0.95;
eps = 0.05;

n = log(1-p)/log(1-1/55/55/pi*eps^2)

%% 7 poisk
p = 0.9;
eps = 0.05;

sz = round(log(1-p)/log(1-1/55/55/pi*eps^2))

%sz = 10^6;
f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);

newsz = round(sz*4/pi);

xy = (rand(2,newsz)-(1/2))*2;
effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
xy = xy(:,effective);

fs = f(xy(1,:),xy(2,:));
[m,i] = min(fs);
x=xy(1,i(1))
y=xy(2,i(1))
%x^2+y^2
mini = m
