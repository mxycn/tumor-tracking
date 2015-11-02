% R =double(imread('50CT.bmp'));
% F =double(imread ('558.bmp'));
% 
% direction=[1 0 0];
% X=[1 1 1];
function [wc3,lambda]=brent(r3d,jx,B,curvature0,stdx)
a=-0.01;
b=0.01;
Epsilon=10^-3;
cgold=0.381966;
% k=1;

IterTimes=500;
% bx=0.5*(a+b);
% if a>b
%     temp=a;
%     a=b;
%     b=temp;
% end

v=a+cgold*(b-a);
w=v;
x=v;
e=0.0;

fx=abs(curvature(r3d,jx)-curvature0);
fv=fx;
fw=fx;
    current_wc3(1)=r3d(jx,2);
    current_wc3(2)=r3d(jx,1);
    current_wc3(3)=stdx-r3d(jx,3);
for iter=1:IterTimes
    xm=0.5*(a+b);
    if abs(x-xm)<Epsilon*2-0.5*(b-a)
      break
    end
    tol1=Epsilon*abs(x)+eps;% tol1=Epsilon*abs(x);
    tol2=2*tol1;
    if abs(e)>tol1    %%%parabola
        r=(x-w)*(fx-fv);
        q=(x-v)*(fx-fw);
        p=(x-v)*q-(x-w)*r;
        q=2*(q-r);
        if q>0
            p=-p;
        end
        q=abs(q);
        etemp=e;
        e=d;
        if not(abs(p)>=abs(0.5*q*etemp)|| p<=q*(a-x)||p>=q*(b-x))
            d=p/q;
            u=x+d;
            if u-a< tol2 || b-u < tol2
                
                d=MySign(Epsilon,xm-x);
            end
        else
            if x>=xm
                e=a-x;
            else
                e=b-x;
            end
            d=cgold*e;
        end
    else
        if x>=xm
               e=a-x;
            else
                e=b-x;
        end
            d=cgold*e;
     end
        
  if abs(d)>=Epsilon
      u=x+d
  else
      u=x+MySign(Epsilon,d)
  end
%     r3d(jx,1:3)
%     wc3(1)=r3d(jx,2);
%     wc3(2)=r3d(jx,1);
%     wc3(3)=stdx-r3d(jx,3);

    wc3=current_wc3+u*B';
    r3d(jx,1)= wc3(2);
    r3d(jx,2)= wc3(1);
    r3d(jx,3)= -wc3(3)+stdx;%将坐标原点移至旋转框架的中心，即世界坐标系的原点
    curvature0
    
  fu=abs(curvature(r3d,jx)-curvature0)
  pause
  if fu<=fx
      if u>=x
          a=x;
      else
          b=x;
      end
      v=w;
      fv=fw;
      w=x;
      fw=fx;
      x=u;
      fx=fu;
  else
      if u<x
          a=u;
      else
          b=u;
      end
      if fu<=fw || w==x
          v=w;
          fv=fw;
          w=u;
          fw=fu;
      else
          if fu<=fv || v==x || v==w
              v=u;
              fv=fu;
          end
      end
  end
end
% minX=x;
% minFx=fx;
% landa=x;
% if abs(x)<=0.001
%     
%     if x<0 
%     
%        x=-0.001;
%     else if x>0 
%        x=0.001;
%     end
%     end
% end
lambda=x;%round(x*1000)/1000;

wc3=r3d(jx,1:3)+x*B';
% Y(1)=round(Y(1));
% Y(2)=round(Y(2));
% Y(3)=round(Y(3)*100)/100;

% Y(1)=round((X(1)+direction(1)*landa)*100)/100;
% Y(2)=round((X(2)+direction(2)*landa)*100)/100;
% Y(3)=round((X(3)+direction(3)*landa)*100)/100;
% fY=fx;
% Y(1)=round((X(1)+direction(1)*x)*10)/10;
% Y(2)=round((X(2)+direction(2)*x)*10)/10;
% Y(3)=round((X(3)+direction(3)*x)*100)/100;

% fY=round(fx*10000)/10000;






