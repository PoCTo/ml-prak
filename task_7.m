%% find root near 
x0 = [-0.357352,0.933969];
cons = @(x)[x(1)^2+x(2)^2-1, []];

[x,fval]=fmincon('objfun',x0,[],[],[],[],[],[],'constraint',...
    optimset('Algorithm','interior-point','TolFun',10^-15,'TolX',10^-17))

%% eps v kruge
real_val = -1.288489227602165;

sz = 10^6;
step = 10^4;
begin = 10^4;



f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);
mins = [];
for s = begin:step:sz
    s
    tic
    newsz = round(sqrt(2)*s);
    xy = ((rand(2,newsz)-(1/2))*2);
    effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
    xy = xy(:,effective);
    xy=xy(:,1:s);


    fs = f(xy(1,:),xy(2,:));
    mn = min(fs)
    mins = [mins, mn];
    toc
end

%%
begin = 10^4;
real_val = -1.288489227602165;

epss = mins - real_val;
%plot(begin:step:sz,epss,'LineWidth',2,'markersize',2)

plot(log(begin:step:sz),-log(epss),'o','markersize',2,'LineWidth',2)
coeffs = polyfit(log(begin:step:sz),-log(epss),1)
hold on
plot(log(begin:step:sz),coeffs(1)*log(begin:step:sz)+coeffs(2),'r','LineWidth',2)
plot(log(begin:step:sz),...
    log( sqrt( (1-(1-0.9).^(log(begin:step:sz).^-1) )/0.000105) ));
set(gca, 'FontSize', 16); 
xlabel('$\mbox{ln}N$','interpreter', 'latex');
ylabel('$-\mbox{ln}\varepsilon$','interpreter', 'latex');
legend('Empirical','Least squares','Location','NorthWest')

% loglog(begin:sz,mins(begin:sz))

%% eps v okr
real_val = -1.288489227602165;

sz = 10^6;
step = 10^4;
begin = 10^4;



f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);
mins = [];
for s = begin:step:sz
    s
    tic
    newsz = round(sqrt(2)*s);
    %xy = ((rand(2,newsz)-(1/2))*2);
    r = rand(1,s)*2*pi;
    xy = [cos(r);sin(r)];
    %effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
    %xy = xy(:,effective);
    %xy=xy(:,1:s);


    fs = f(xy(1,:),xy(2,:));
    mn = min(fs)
    mins = [mins, mn];
    toc
end

%%
begin = 10^4;
real_val = -1.288489227602165;

epss = mins - real_val;
%plot(begin:step:sz,epss,'LineWidth',2,'markersize',2)

plot(log(begin:step:sz),-log(epss),'o','markersize',2,'LineWidth',2)
coeffs = polyfit(log(begin:step:sz),-log(epss),1)
hold on
plot(log(begin:step:sz),coeffs(1)*log(begin:step:sz)+coeffs(2),'r','LineWidth',2)

set(gca, 'FontSize', 16); 
xlabel('$\mbox{ln}N$','interpreter', 'latex');
ylabel('$-\mbox{ln}\varepsilon$','interpreter', 'latex');
legend('Empirical','Least squares','Location','NorthWest')

% loglog(begin:sz,mins(begin:sz))