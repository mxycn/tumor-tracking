function [bb rr]=six_var_regress(n,x1,x2,x3,y1,y2,y3,z1,z2,z3)
% each z = b1 + b2*x1 + b3*x2 + b4*x3 + b5*y1 +  b6*y2 + b7*y3
[bb1,bint1,r1,rint1,stats1] = regress(z1,[ones(n,1) x1 x2 x3 y1 y2 y3]);
[bb2,bint2,r2,rint2,stats2] = regress(z2,[ones(n,1) x1 x2 x3 y1 y2 y3]);
[bb3,bint3,r3,rint3,stats3] = regress(z3,[ones(n,1) x1 x2 x3 y1 y2 y3]);
rr = max([stats1 stats2 stats3]);
bb = [bb1 bb2 bb3];
