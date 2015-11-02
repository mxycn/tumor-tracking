function wc3=back_proj_new(centerpoint,a,b,c,y_min,y_max,r,f,gAng2,u,v)
%
% % see Appendices of Liuwu for details
% Programmed by: Huagang Yan
% 07/15/2010
%% Input
% k: the slope of the tracks in x-y plane in the world coordinates,
% obtained in the first 4 seconds by regression based on the triangulation
% between kV data and MV data.
% r: the source to axis distance(SAD)
% f: the source to imager plate distance
% gAng1: the gantry angle of the last projection
% gAng2: the gantry angle of the present projection
% u1,u2: the pixel coorinate x on the portal imager 
%        left + (data coordinate y, liuwu coordinate x)
% v1,v2: the pixel coorinate z on the portal imager  
%        infer + (data coordinate x, liuwu coordiante y)
%   Anterior+: data coordinate z+, Liuwu coordinate z-
%--------------------------------------------------------------------------

O=[-r*sin(gAng2*pi/180.0),0,r*(1-cos(gAng2*pi/180.0))]';
Ap=[u,v,f]';
A=[0 0 0]';
R=[cos(gAng2*pi/180.0) 0 sin(gAng2*pi/180.0); 0 1 0; -sin(gAng2*pi/180.0) 0 cos(gAng2*pi/180.0)];

A=R*Ap+O; % the 3D coordinate of the current projection point

B=A-O;
is_useyz=0;

% in the following, we are going to seek a point on the current projection
% line that is of the shortest distance extended from the the last
% projection line in the plane that satisfies y= k x.-- according to liuwu's convention
% intermediate step1
% in liuwu's coordinate, x=k1(1)*y+k1(2);z=k2(1)*y+k2(2);z=k3(1)*x+k3(2);
% and the new form is:
% (x-centerpoint(2))/b=(y-centerpoint(1))/a=(z-std+centerpoint(3))/-c

k1(1)=b/a;
k1(2)=-b/a*centerpoint(1)+centerpoint(2);
k2(1)=-c/a;
k2(2)=c/a*centerpoint(1)+r-centerpoint(3);

 a1=-k2(1)*(-k1(2)+O(1))+k1(1)*(-k2(2)+O(3));
 b1=-k2(1)*B(1)+k1(1)*B(3);
 a2=-k1(2)+O(1)-k1(1)*O(2);
 b2=B(1)-k1(1)*B(2);
 a3=k2(2)-O(3)+k2(1)*O(2);
 b3=-B(3)+k2(1)*B(2);
 
 s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
%  if (abs(k1(1))>2.0 & abs(k2(1))>2.0) | (abs(k1(1))>3.0 & abs(k2(1))>1.0) | (abs(k1(1))>1.0 & abs(k2(1))>3.0)
%      s=(O(3)-k3(1)*O(1)-k3(2))/(B(1)*k3(1)-B(3));
%      is_useyz=1;
%  end
     
P=[0 0 0]';
% the 3D coordinates of the current target.
P=O+s*B;
margin=0.0;
% marginx=1.0;
% if is_useyz==0 
    % if the point is close to the check line extension, instead of the
    % checkline segment, then use another method.
    bb=P(2)-k1(1)*(k1(2)-P(1))-k2(1)*(k2(2)-P(3));
    aa=k1(1)*k1(1)+k2(1)*k2(1)+1;    
    y=bb/aa; % the y-coordinate of the point on the trajectory line that is closest to the target 
 % negative means conservative, positive means following trend     
%     if y_min-y>margin     
    if y_min-y>margin*abs(a)             
        xmin=k1(1)*y_min+k1(2);
        ymin=y_min;
        zmin=k2(1)*y_min+k2(2);
        a1=O(1)-xmin;
        a2=O(2)-ymin;
        a3=O(3)-zmin;
        b1=B(1);
        b2=B(2);
        b3=B(3);
        s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
        P=O+s*B;
    elseif y-y_max>margin*abs(a)
%     elseif y-y_max>margin 
        xmax=k1(1)*y_max+k1(2);
        ymax=y_max;
        zmax=k2(1)*y_max+k2(2);
        a1=O(1)-xmax;
        a2=O(2)-ymax;
        a3=O(3)-zmax;
        b1=B(1);
        b2=B(2);
        b3=B(3);
        s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
        P=O+s*B;    
    end    
    err=sqrt((a1+b1*s)^2+(a2+b2*s)^2+(a3+b3*s)^2)/sqrt(k1(1)*k1(1)+k2(1)*k2(1)+1);
wc3(1:3)=P;
wc3(4)=err;
return;