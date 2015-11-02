function WaperR_adjusted=periodic_evaluate(s,T,f)
l=length(s);
m=round(T*f);
n=floor(l/f/T);
A=zeros(n,m);
for i=0:n-1
    A(i+1,:)=s(round(i*T*f)+1:round(i*T*f)+m);
end
%% SVD
% AA=A';
% xx=1:length(A(:));
% plot(xx,AA(:));
% pause
[u,sigma,v]=svd(A);
e=diag(sigma);
Waper=(norm(e).^2-e(1).^2)/norm(e).^2;
p=sum(e~=0);
WaperR=p/(p-1)*Waper;
alpha=u(:,1)*e(1)*v(:,1)';
EE=sum(alpha.^2,2);
ee=EE/sum(EE);
ehr=-sum(ee.*log(ee)/log(n));
WaperR_adjusted=1-ehr*(1-WaperR)
%%2015-0801 by Mxy
