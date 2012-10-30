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

%% 1.2 ��� for Bi
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

%% 1.2 ��� for nBi
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

%% 1.3 �������
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

%% 4.1 cauchy
x0=0;
ggamma = 2;
sz = 5000;
fontsize = 16;
bars = 100;
leftcons = -10;
rightcons = 10;

seq = arrayfun(@(x)generateCauchyRandomValue(x0, ggamma),1:sz);
seq = seq(find(seq>=(leftcons+x0) & seq<=(rightcons+x0)));
[f,n] = hist(seq,bars);
bar(n,f/trapz(n,f));
hold on
h=ezplot(@(x)1/pi/ggamma/(1+(x-x0)^2/ggamma^2),[min(seq),max(seq)]);
set(h,'LineWidth',2, 'Color', 'r')
set(gca,'YLim',[0,1/pi/ggamma/(1+(x0-x0)^2/ggamma^2)]*1.1);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% 4.2 Normal
sigma2 = 1;
mu = 0;
sz = 5000;
fontsize = 16;
bars = 100;

seq = arrayfun(@(x)generateNormalRandomValueVonNeuman(mu,sigma2),1:sz);
seq = seq(find(seq ~= Inf));
[f,n] = hist(seq,bars);
bar(n,f/trapz(n,f));
hold on
h=ezplot(@(x)1/sqrt(sigma2*2*pi)*exp(-(x-mu)^2/2/sigma2),[min(seq),max(seq)]);
set(h,'LineWidth',2, 'Color', 'r')
%set(gca,'YLim',[0,l]);
title '';
set(gca, 'FontSize', fontsize);
xlabel '';

%% 4.3 speeds
xs=[500:10:5000];
times1=[];
times2=[];
for i=xs
    i
    tic;
    cur = 0;
    while (cur<i)
        if (generateNormalRandomValueVonNeuman(0,1)~=Inf)
            cur = cur+1;
        end
    end
    times1 = [times1, toc];
    
    tic;
    cur = 0;
    while (cur<i)
        t=[generateNormalRandomValuePair1(0,1),...
            generateNormalRandomValuePair2(0,1)];
        cur = cur+2;
        
    end
    times2 = [times2, toc];
end

%% 4.3 speedgraph
fontsize = 16;

plot(xs(5:end),times1(5:end),'r','LineWidth',2);
hold on;
plot(xs(5:end),times2(5:end),'g','LineWidth',2);
xlabel('$N$','interpreter', 'latex');
ylabel('$T(N)$','interpreter', 'latex');
h=legend('Von Neuman','Polar');
set(h,'interpreter','latex','Location','NorthWest');
set(gca, 'FontSize', fontsize);
xlim([xs(5),xs(end)]);
%% 4.3 otnoshenie
fontsize = 16

plot(xs(5:end),times1(5:end)./times2(5:end),'r','LineWidth',2);
xlabel('$N$','interpreter', 'latex');
%ylabel('$T_{VN}/T_{Polar}$','interpreter', 'latex');
h=legend('$T_{Neuman}/T_{Polar}$','Polar');
set(h,'interpreter','latex','Location','SouthWest');
set(gca, 'FontSize', fontsize);
xlim([xs(5),xs(end)]);

%% 5.1 zbch-normal
sigma2 = 1;
mu = 0;
sz = 10000;
fontsize = 16;
bars = 100;

seq1 = arrayfun(@(x)generateNormalRandomValuePair1(mu,sigma2),1:sz);
seq2 =arrayfun(@(x)generateNormalRandomValuePair2(mu,sigma2),1:sz);

sum1=cumsum(seq1);
sum2=cumsum(seq2);
sum=(sum1+sum2)./[2:2:(sz*2)];
plot(80:2:(sz*2),sum(40:end),'b','LineWidth',1.5);
hold on;
plot([80,sz*2],[mu,mu],'r','LineWidth', 2);
xlabel('$n$','interpreter', 'latex');
ylabel('$S_n$','interpreter', 'latex');
set(gca, 'FontSize', fontsize);

%% 5.2 no-zbch-cauchy
p_a = 0;
p_b = 1;
sz = 100000;
fontsize = 16;
bars = 100;

seq = p_a+p_b*arrayfun(@(x)generateCauchyRandomValue(0,1),1:sz);

sum=(cumsum(seq))./[1:(sz)];
plot(40:(sz),sum(40:end),'b','LineWidth',1.5);
hold on;
plot([40,sz],[p_a,p_a],'r','LineWidth', 2);
xlabel('$N$','interpreter', 'latex');
ylabel('$S_n$','interpreter', 'latex');
set(gca, 'FontSize', fontsize); 
xlim([40,sz]);

%% save
saveas(gca, 'D:\Sync\Dropbox\ml\img\file1.eps', 'psc2')

%% sandbox
generateCauchyRandomValue(1,1)