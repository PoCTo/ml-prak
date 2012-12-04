%% 7 (eps,p) -> n
p = 0.95;
eps = 0.05;

n = log(1-p)/log(1-1/55/55/pi*eps^2)

%% 7 poisk na kruge
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
x^2+y^2
mini = m

%% eps graph
n = 10^6;
f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);

sz=n;
newsz = round(sz*4/pi);

xy = (rand(2,newsz)-(1/2))*2;
effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
xy = xy(:,effective);

fs = f(xy(1,:),xy(2,:));
[m,i] = min(fs);
x=xy(1,i(1))
y=xy(2,i(1))
x^2+y^2
mini = m


mins=[];
mins=[0]
for i=10^3:10^3:sz
    i
    mins=[mins,min(mins(end),min(fs(i-10^3+1:i)))];
end

%% plott
empeps = abs(mins+1.28848922760216498961553789);
plot(10^3:10^3:sz,empeps(2:end),'r','LineWidth',2);
hold on
ns=10^3:10^3:sz;
p=0.9;
plot(10^3:10^3:sz,sqrt(ns.^-1)*sqrt(-log(1-p)/0.000105),'LineWidth',2)

legend('Empirical','Upper Bound')
fontsize = 16;
set(gca, 'FontSize', fontsize); 


%% 7 orkuzhnost
p = 0.95;
eps = 0.05;

n = log(1-p)/log(1-1/55/pi*eps)

%% 7 poisk na okr
p = 0.999;
eps = 0.001;

sz = round(log(1-p)/log(1-1/55/pi*eps))

f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);

r = rand(1,sz)*2*pi;
xy = [cos(r);sin(r)];
%effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
%xy = xy(:,effective);

fs = f(xy(1,:),xy(2,:));
[m,i] = min(fs);
x=xy(1,i(1))
y=xy(2,i(1))
x^2+y^2
mini = m

%% eps graph okr
n = 10^6;
f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);

sz = n

f=@(x,y)(x.^3).*sin(x.^-1)+10*x.*(y.^4).*cos(y.^-1);

r = rand(1,sz)*2*pi;
xy = [cos(r);sin(r)];
%effective = find(xy(1,:).^2+xy(2,:).^2 <= 1);
%xy = xy(:,effective);

fs = f(xy(1,:),xy(2,:));
[m,i] = min(fs);
x=xy(1,i(1))
y=xy(2,i(1))
x^2+y^2
mini = m


mins=[];
mins=[0]
for i=10^3:10^3:sz
    i
    mins=[mins,min(mins(end),min(fs(i-10^3+1:i)))];
end

%% plott
empeps = abs(mins+1.28848922760216498961553789);
plot(10^3:10^3:sz,empeps(2:end),'r','LineWidth',2);
hold on
ns=10^3:10^3:sz;
p=0.9;
plot(10^3:10^3:sz,(ns.^-1)*(-log(1-p)/0.0057),'LineWidth',2)

legend('Empirical','Upper Bound')
fontsize = 16;
set(gca, 'FontSize', fontsize); 



%% 8 eps |-> N
eps = 0.05;
f = @(x,y)x.^2;

N = round(9*1/eps^2)

%% 8 result
rng('shuffle');

grid_sz = 30;
experiment_count = 100;

sz=0;

f = @(x,y)x^2;

[xs,ys] = meshgrid(-1:(1/grid_sz):1);
fs = Inf(size(xs));

i=1;
maxj=0;

for x = 1:size(xs,1)
    x
    for y = 1:size(xs,2)
        if ((xs(x,y))^2+(ys(x,y))^2>1) 
            continue;
        end
        sum=0;
        for i=1:experiment_count
            x1 = x; y1=y;
            %'begin'
            j = 0;
            while ((xs(x1,y1))^2+(ys(x1,y1))^2 < 1) 
                r=randi(4);                
                if (r == 1)
                    x1 = x1+1;
                elseif (r==2)
                    x1 = x1-1;
                elseif (r==3)
                    y1 = y1+1;
                else
                    y1 = y1-1;
                end                    
                %[x1] = [x1] +[r-1]*(randi(2)/2-1/2)*2/grid_sz;
                %[y1] = [y1] +[2-r]*(randi(2)/2-1/2)*2/grid_sz;                              
            %    [x1,y1]   
                j=j+1;
            end
            if (j>maxj)
                maxj=j;
            end
            %'end'
            %[xs(x1,y1), ys(x1,y1), (xs(x1,y1))^2+(ys(x1,y1))^2]
            norm = (xs(x1,y1)^2+ys(x1,y1)^2);
            sum=sum+f(xs(x1,y1)/norm,ys(x1,y1)/norm);
            %xs(x1,y1)
        end
        %[sum, xs(x,y), ys(x,y)]
        fs(x,y) = sum/experiment_count;
        %xs(i) = x;
        %ys(i) = y;
    end
    
end

maxj

%% 8 - surf
surf(xs,ys,fs)

%% 8 - real
f_s = @(x,y)(1/2)*(x.^2-y.^2+1);
surf(xs,ys,f_s(xs,ys));

%% 9 - Wiener
set(0,'RecursionLimit',5000)

eps = 0.00001;
a = 0.5; 

val1 = generateNormalRandomValuePair1(0,1,rand(1),rand(1));

[X,Y] = findWienerRec(0, 1,...
    0, val1,...
    eps, a);
%length(find(X(1:end-1)-X(2:end)>eps))
%% 9 - plot
fontsize = 16;

plot(X,Y,'o','markersize',1)
hold on
set(gca, 'FontSize', fontsize); 
xlabel('$t$','interpreter', 'latex');
ylabel('$W_t$','interpreter', 'latex');

%% 10 - Ulenbeck
set(0,'RecursionLimit',5000)

eps = 0.0001;
a = 0.8; 
sigma2 = 0.1;
lambda = 10^-5;

val0 = generateNormalRandomValue(0,1,rand(1),rand(1))
val1 = generateNormalRandomValue(exp(-lambda)*val0,...
    sigma2-sigma2*exp(-2*lambda),rand(1),rand(1))

[X,Y] = findOrnsteinRec(0, 1,...
    val0, val1,...
    eps, a, lambda, sigma2)
%length(find(X(1:end-1)-X(2:end)>eps))
%% 10 - Plot
fontsize = 16;

plot(X,Y,'o','markersize',1)
hold on
set(gca, 'FontSize', fontsize); 
xlabel('$t$','interpreter', 'latex');
ylabel('$X(t)$','interpreter', 'latex');

%% 11.1 SMO
fontsize = 16;
lambda = 0.2;
Tmax = 500;

t=0; ts=[0];
while(t<Tmax)
    t=t+generateExponentialRandomValue(lambda);
    ts=[ts, t];
end
ts=ts(1:end-1);

xs=zeros(1,length(ts));
for i=1:length(ts)
    xs(i)=generateChi2RandomValue(10);
end

ends=zeros(1,length(ts));
ends(1)=xs(1)+ts(1);
for i=2:length(ts)
    ends(i) = xs(i) + max(ts(i),ends(i-1));
end
all = [ts,ends;ones(1,length(xs)),-ones(1,length(xs))];
sorted = sortrows(all');
%sorted(:,2)=cumsum(sorted(:,2));
sorted = sorted(find(sorted(:,1)<=Tmax),:)

plot([ts(1),ts(1)],[0,1])
k=1;
hold on;
for i=2:length(sorted(:,1))
    plot([sorted(i-1,1),sorted(i,1)],[k,k]);
    plot([sorted(i,1),sorted(i,1)],[k,k+sorted(i,2)]);
    k=k+sorted(i,2);
end
set(gca, 'FontSize', fontsize); 
ylabel('Queue Size','interpreter', 'latex');
xlabel('$t$','interpreter', 'latex');

%% 11.2 - Insurance
fontsize = 16;
lambda = 0.5;
k = 0.5;
x_m=1;
u0 = 10;
c=-1;
Tmax = 500;

t=0; ts=[0];
while(t<Tmax)
    t=t+generateExponentialRandomValue(lambda);
    ts=[ts, t];
end
ts=ts(1:end-1);

%xs=zeros(1,length(ts));
xs=generateParetoRandomMatrix(x_m,k,1,length(ts))

u=u0;
plot([0],u0);
hold on;
for i=2:length(ts)
    plot([ts(i-1),ts(i)],[u,u + c*(ts(i)-ts(i-1))]);
    if (u + c*(ts(i)-ts(i-1)) - xs(i) < 0)
        plot([ts(i),ts(i)],[u + c*(ts(i)-ts(i-1)),0]);
        break
    else
        plot([ts(i),ts(i)],...
            [u + c*(ts(i)-ts(i-1)),u + c*(ts(i)-ts(i-1)) - xs(i)]);
    end
    u = u + c*(ts(i)-ts(i-1)) - xs(i);    
end
set(gca, 'FontSize', fontsize); 
xlabel('$u(t)$','interpreter', 'latex');
ylabel('$t$','interpreter', 'latex');
%% 11.2 - razorenie
fontsize = 16;
lambda = 0.5;
k = 0.5;
x_m=1;
u0 = 10;
c=-1;
Tmax = 500;

bars = 40;
sz = 1000;
tic
ruins = arrayfun(@(x)razorenie(lambda,u0,c,k,x_m,Tmax),1:sz);
toc

leninf=length(find(ruins==Inf));
percinf=leninf/length(ruins)
ruins = ruins(find(ruins~=Inf));

[f,n] = hist(ruins,bars);
%bar(n,f/trapz(n,f))
hist(ruins,bars)
%% save
saveas(gca, 'D:\Sync\Dropbox\ml\img\file1.eps', 'psc2')

%% save linux
saveas(gca, '/home/pocto/ml-new/ml-prak/tex/eps/loglog_krug.eps','psc2')

