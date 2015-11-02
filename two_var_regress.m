function [bb rr]=two_var_regress(n,x,y,z)
% z = a + b*x + c*y
[bb,bint,r,rint,stats] = regress(z,[ones(n,1) x y] );
rr = stats(1);
return