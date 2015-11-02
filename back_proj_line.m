function [wc3 B cost]=back_proj_line(centerpoint,kp,prev,prep,marginx,y_min,y_max,r,f,gAng2,u,v)
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

k1(1)=kp(2)/kp(1);
k1(2)=-kp(2)/kp(1)*centerpoint(1)+centerpoint(2);
k2(1)=-kp(3)/kp(1);
k2(2)=kp(3)/kp(1)*centerpoint(1)+r-centerpoint(3);

%  a1=-k2(1)*(-k1(2)+O(1))+k1(1)*(-k2(2)+O(3));
%  b1=-k2(1)*B(1)+k1(1)*B(3);
%  a2=-k1(2)+O(1)-k1(1)*O(2);
%  b2=B(1)-k1(1)*B(2);
%  a3=k2(2)-O(3)+k2(1)*O(2);
%  b3=-B(3)+k2(1)*B(2);
%  
%  s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
P=[0 0 0]';
% the 3D coordinates of the current target.
factorv = 0.02;
[s cost] = find_s(O,B,centerpoint,kp,prev,prep,marginx,y_min,y_max,r,factorv);
P=O+s*B;
% if is_useyz==0 
% err=sqrt((a1+b1*s)^2+(a2+b2*s)^2+(a3+b3*s)^2)/sqrt(k1(1)*k1(1)+k2(1)*k2(1)+1);
wc3(1:3)=P;
% wc3(4)=err;
return;