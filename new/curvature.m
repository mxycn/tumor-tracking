function k=curvature(s,j)
if  j>=3
    v=s(j,:)-s(j-1,:);
    a=s(j,:)-2*s(j-1,:)+s(j-2,:);
    k=sqrt(norm(a).^2.*norm(v).^2-(v*a').^2);
else k=0;
end