pogreshn = 0
n = 4
x1 = cos((2*1-1)*pi/2/n)
x2 = cos((2*2-1)*pi/2/n)
x3 = cos((2*3-1)*pi/2/n)
x4 = cos((2*4-1)*pi/2/n)

a=x1
b=x2
c=x3
d=x4

xs=linspace(-0.999,0.999,200);
%ys=arrayfun(@(x)h(x),xs)

%h = @(x)1/sqrt(1-x^2)*(x-x1)*(x-x2)*(x-x3)*(x-x4)
int = ((3 + 4* c* d + 4* b* (c + d) + 4*a* (b + c + d + 2* b* c *d))*pi)/8