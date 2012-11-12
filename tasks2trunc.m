%% 8 eps |-> N
eps = 0.05;
f = @(x,y)x.^2;

N = round(9*1/eps^2)

%% 8 result

grid_sz = 20;
experiment_count = 20;

sz=0;


for x = -1:(1/grid_sz):1
    x
    sz=sz+length((x^2-1):(1/grid_sz):(1-x^2));
end

fs=zeros(1,sz);
xs=zeros(1,sz);
ys=zeros(1,sz);

i=1;
for x = -1:(1/grid_sz):1
    x
    for y = (x^2-1):(1/grid_sz):(1-x^2)
        for i=1:experiment_count
            x1 = x; y1=y;
            sum=0;
            while (x1^2+y1^2<1) 
                r=randi(2);
                [x1] = [x1] +[r-1]*(randi(2)/2-1/2)*2/grid_sz;
                [y1] = [y1] +[2-r]*(randi(2)/2-1/2)*2/grid_sz;
            
            end
            sum=sum+f(x1,y1);
        end
        fs(i) = sum/experiment_count;
        xs(i) = x;
        ys(i) = y;
    end
    
end

%%
length(-1:(1/grid_sz):1)
length((x^2-1):(1/grid_sz):(1-x^2))
length(fs)
reshape(fs, length(-1:(1/grid_sz):1),length((x^2-1):(1/grid_sz):(1-x^2)))

surf(xs,ys,fs)