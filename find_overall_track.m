function [linearity,range_o3d]=find_overall_track(j0,allpoints,ko3d)
% compared to initializing, this function adopts PCA method.

% to initialize the tracking by turning on kV imager. I.e.,
% to use the fisrt 'runlength' true points to obtain the early 
% correlation and the checkpoint position, plus the 
% first 'runlength' reconstructed positions.
% the series starts from j0, 

% linearity: an index to indicate linearity of the trajectory.
% a,b,c:   When is_line=0. They are used jointly with n_domin,
%          if n_domin=1, the fit plane is x=a*y+b*z+c.
%          When is_line=1. They are used jointly with centerpoint, the most
%          fit line is:
%          (x-centerpoint(1))/a=(y-centerpoint(2))/b=(z-centerpoint(3))/c
% All data come from the measured data r3d.

o3dx=zeros(allpoints,3);
range_o3d=[0 0 0];
for jx=1:allpoints
    %% adding error to the reconstruction(estimation)        
    o3dx(jx,1)=ko3d(jx+j0-1,1);
    o3dx(jx,2)=ko3d(jx+j0-1,2);
    o3dx(jx,3)=ko3d(jx+j0-1,3);     
end
[coeffo,scoreo,rootso] = princomp(o3dx);
linearity = rootso(1)/sum(rootso);
xx=range(ko3d(j0:j0+allpoints-1,1:3),1);
range_o3d=xx(1,:);
return