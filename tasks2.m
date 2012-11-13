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

%% 2 - surf
surf(xs,ys,fs)

%% 2 - real
f_s = @(x,y)(1/2)*(x.^2-y.^2+1);
surf(xs,ys,f_s(xs,ys));

%% save
saveas(gca, 'D:\Sync\Dropbox\ml\img\file1.eps', 'psc2')

%% save linux
saveas(gca, '/home/pocto/ml-prak/tex/eps/file.eps','psc2')

