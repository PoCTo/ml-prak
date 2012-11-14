function [X,Y]=findWienerRec(x,y,x_val,y_val,eps,a)
    if (abs(y-x)<eps) 
        X=[];
        Y=[];
        return
    end
    %a = rand(1);
    new_val = generateNormalRandomValue((1-a)*x_val+a*y_val,...
        a*(1-a)*(y-x),rand(1),rand(1));
    %new_val = generateNormalRandomValueVonNeumanTrue((1-a)*x_val+a*y_val,...
    %    a*(1-a)*(y-x));
    [Xless,Yless] = findWienerRec(x,x+(y-x)*a,x_val,new_val,eps,a);
    [Xmore,Ymore] = findWienerRec(x+(y-x)*a,y,new_val,y_val,eps,a);
    Y = [Yless,...
         new_val,...
         Ymore];
    X = [Xless, x+(y-x)*a , Xmore];
    
