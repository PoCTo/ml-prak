%% 1.1 bi(n,p)
n = 20;
p = 0.9;
sz = 1000;
fontsize = 16;

hist(generateBinomialRandomSequence(n, p, sz),[0:n]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','none','EdgeColor','black', 'LineWidth',2);
set(gca, 'FontSize', fontsize);
hold on
right = generateRightBinomialCount(n, p, sz)
bar(0:n,right,'FaceColor','none','EdgeColor','b', 'LineWidth',2);
%hist(generateBinomialRandomSequence(100, 0.5, 5000),[0:100])

%% 1.1 nbi(n,p)
r = 10;
p = 0.8;
sz = 1000;
fontsize = 16;

seq = generateNegativeBinomialRandomSequence(r, p, sz);
hist(seq, [min(seq):max(seq)]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','none','EdgeColor','black', 'LineWidth',2);
set(gca, 'FontSize', fontsize);
hold on
right = generateRightNegativeBinomialCount(r, p, sz)
bar(0:length(right)-1,right,'FaceColor','none','EdgeColor','b', 'LineWidth',1.5);
%hist(generateBinomialRandomSequence(100, 0.5, 5000),[0:100])

%% 1.2 ÇÁ× for Bi
sz = 5000;
n = 1;
p = 0.5;
fontsize = 16;

seq = generateBinomialRandomSequence(n, p, sz);
empE = cumsum(seq);
empE = empE(20:end)./[20:sz];
plot([20:sz],empE);
hold on
plot([20:sz],n*p,'Color','r')
set(gca, 'FontSize', fontsize);

%% 1.2 ÇÁ× for nBi
sz = 10000;
r = 4;
p = 0.3;
fontsize = 16;

seq = generateNegativeBinomialRandomSequence(r, p, sz);
empE = cumsum(seq);
empE = empE(50:end)./[50:sz];
plot([50:sz],empE);
hold on
plot([50:sz],r*(1-p)/p,'Color','r')
set(gca, 'FontSize', fontsize);

%% 1.3 Îðëÿíêà
N = 1000
p = 0.5
fontsize = 16;

seq = generateBernulliRandomSequence(p, N+1);
seq(find(seq==0))=-ones(1,length(find(seq==0)));
xs = (0:N)/N;
ys = cumsum(seq)/sqrt(N);
plot(xs,ys);
set(gca, 'FontSize', fontsize);

%% 2.4
prec = 50;
sz = 5000;
fontsize = 16;

seq = generateKantorRandomSequence(prec, sz);
plotEmpiricalDistribution(seq,'b');
set(gca, 'FontSize', fontsize);
xlabel('$x$', 'interpreter', 'latex');
ylabel('$F(x)$', 'interpreter', 'latex');

%% 2.3 symm
prec = 50;
sz = 5000;
fontsize = 16;

seq = generateKantorRandomSequence(prec, sz);
plotEmpiricalDistribution(seq, 'b');
hold on;
plotEmpiricalDistribution(1-seq, 'r--');
set(gca, 'FontSize', fontsize);
xlabel('$f(x)$', 'interpreter', 'latex');
ylabel('$F(f(x))$', 'interpreter', 'latex');
h=legend('$f(X)=X$','$f(X)=1-X$');
set(h,'interpreter','latex','Location','NorthWest');

%% 2.3 self-sim
prec = 50;
sz = 5000;
fontsize = 16;

seq = generateKantorRandomSequence(prec, sz);
plotEmpiricalDistribution(seq, 'b');
hold on;
plotEmpiricalDistribution(1-seq, 'r--');
set(gca, 'FontSize', fontsize);
xlabel('$f(x)$', 'interpreter', 'latex');
ylabel('$F(f(x))$', 'interpreter', 'latex');
h=legend('$f(X)=X$','$f(X)=(X| X<\frac 13)$');
set(h,'interpreter','latex','Location','NorthWest');

%% 3.1 exp
l = 5;
sz = 1000;
fontsize = 16;
bars = 100;

seq = arrayfun(@(x)generateExponentialRandomValue(l),1:sz);
[f,n] = hist(seq,bars)
bar(n,f/trapz(n,f));
hold on
h=ezplot(@(x)(l*exp(-l*x)),[min(seq),max(seq)]);
set(h,'LineWidth',2, 'Color', 'r')
set(gca,'YLim',[0,l]);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% 3.1 poiss
l = 10;
sz = 1000;
fontsize = 16;
bars = 100;

seq = arrayfun(@(x)generatePoissonRandomValue(l),1:sz);
hist(seq,min(seq):max(seq))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black', 'LineWidth',2);
hold on
generateRightPoissonCount(l,sz)
bar(0:(2*l+2-1),generateRightPoissonCount(l,sz),'FaceColor','none','EdgeColor','b', 'LineWidth',1.5)

%set(h,'LineWidth',2, 'Color', 'r')
%set(gca,'YLim',[0,l]);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% 3.2 Normal
sigma2 = 5;
mu = -5;
sz = 5000;
fontsize = 16;
bars = 100;
seq = arrayfun(@(x)generateNormalRandomValuePair1(mu,sigma2),1:sz);
seq = [seq,arrayfun(@(x)generateNormalRandomValuePair2(mu,sigma2),1:sz)];
[f,n] = hist(seq,bars)
bar(n,f/trapz(n,f));
hold on
h=ezplot(@(x)1/sqrt(sigma2*2*pi)*exp(-(x-mu)^2/2/sigma2),[min(seq),max(seq)]);
set(h,'LineWidth',2, 'Color', 'r')
%set(gca,'YLim',[0,l]);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% 3.3 chi2
k=5;
sz = 5000;
fontsize = 16;
bars = 100;

seq = arrayfun(@(x)generateChi2RandomValuePair1(k),1:sz);
seq = [seq,arrayfun(@(x)generateChi2RandomValuePair2(k),1:sz)];
[f,n] = hist(seq,bars)
bar(n,f/trapz(n,f));
hold on
cnst=(1/2)^(k/2)/gamma(k/2);
h=ezplot(@(x)cnst*x^(k/2-1)*exp(-x/2),[min(seq),max(seq)]);
set(h,'LineWidth',2, 'Color', 'r')
%set(gca,'YLim',[0,l]);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% save
saveas(gca, 'D:\dev\ml-prak\img\file1.eps', 'psc2')

%% sandbox
generateRightPoissonCount(10,sz)
