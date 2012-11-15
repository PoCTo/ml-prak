function t=razorenie(lambda,u0,c,k,x_m,Tmax)

%lambda = 0.5;
%k = 0.2;
%x_m=0.2;
%u0 = 40;
%c=12;
%Tmax = 500;

t=0; ts=[0];
while(t<Tmax)
    t=t+generateExponentialRandomValue(lambda);
    ts=[ts, t];
end
ts=ts(1:end-1);

%xs=zeros(1,length(ts));
xs=generateParetoRandomMatrix(x_m,k,1,length(ts));

u=u0;
for i=2:length(ts)
    if (u+c*(ts(i)-ts(i-1))<0)
        t=ts(i-1)+u/abs(c);
        return
    end
    if (u + c*(ts(i)-ts(i-1)) - xs(i) < 0)
        t=ts(i);
        return
    end
    u = u + c*(ts(i)-ts(i-1)) - xs(i);    
end
t=Inf;