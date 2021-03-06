function [r3dx,is_discon,centerx,kx,ko,range_o3d,is_jump,linearity]=proj_backproj_fixed(seg_err,...
    ko3d,r3d,jx,j0,gAng,runlength,sid,stdx);
% this function projects the actual points to the MV imager and 
% reconstruct them

% r3dx returns the new reconstruction result.
% is_discon signals the discontinuity in the original data.
% is_jump signals the large velocity in the actual 3D motion.

% seg_err is the detection error(=0.5mm)
% k is the no. of fraction.
% o3d is the actual 3D positions
% r3d is the estimated 3D positions
% jx is the no. of point to be estimated.
% j0 is the starting number of the current trajectory
% delta_angle is the angle of gantry rotation during the first 4 seconds.

% Initialize the flags. 
is_discon=0;
centerx=[0 0 0];
kx=[0 0 0];
ko=[0 0 0];
ym=[0 0];
is_jump=0;
isplane=0;
linearity=1;
xfit=zeros(jx-j0,3);
r3d0=zeros(jx-j0,3);
r3dd=zeros(jx-j0,3);
ko3dd=zeros(jx-j0,3);
% the current gantry angle.

% if the velocity is too large, abandon the current trial
if abs(ko3d(jx,1)-ko3d(jx-1,1)) > 6.5 | ...
        abs(ko3d(jx,2)-ko3d(jx-1,2)) > 6.5 | ...
        abs(ko3d(jx,3)-ko3d(jx-1,3)) > 6.5
    is_jump=1; 
    is_discon=2;
    r3dx(1)= 0;
    r3dx(2)= 0;
    r3dx(3)= 0;
    range_o3d=[0 0 0];    
    return
end
% it has been compared to 1:jx-j0. The following is better.
for ix=1:runlength     
    r3d0(ix,1)=r3d(j0+ix-1,2);
    r3d0(ix,2)=r3d(j0+ix-1,1);
    r3d0(ix,3)=stdx-r3d(j0+ix-1,3);
end

kki=1;
for kkk=1:jx-j0
    if ~isnan(r3d(kkk+j0-1,1)) 
        r3dd(kki,1:3)=r3d(kkk+j0-1,1:3);
        ko3dd(kki,1:3)=ko3d(kkk+j0-1,1:3);
        kki=kki+1;
    end
end

p2d=[0 0];
wc3=[0 0 0];
% if the point data in the database is continuous(based on the time 
% information), perform estimation, flag it otherwise.
if (ko3d(jx+1,4)-ko3d(jx,4)<0.2)
    o_3d=[0 0 0 0];
    o_3d(:)=ko3d(jx,:);
    p2d(1:2)=proj([o_3d(1) o_3d(2) o_3d(3)],gAng,sid,stdx);
    [coeff,score,roots] = princomp(r3dd(1:kki-1,1:3));
    [coeffo,scoreo,rootso] = princomp(ko3dd(1:kki-1,1:3));
    centerx = mean(r3dd(1:kki-1,1:3),1);
    xfit=repmat(centerx,kki-1,1) + score(:,1)*coeff(:,1)';
    linearity = roots(1)/sum(roots);
    if linearity>0.07
        isplane=0;
    else
        isplane=1;
    end        
    dir_vector = coeff(:,1)';
    diro_vector = coeffo(:,1)';
    kx(:)=dir_vector(:);
    ko(:)=diro_vector(:);
    %% add segmentation error to the data
    p2d(1)=p2d(1) + seg_err * randn();
    p2d(2)=p2d(2) + seg_err * randn();
    ymin=min(xfit(:,1));
    ymax=max(xfit(:,1));  
    ym=[ymin ymax];
    if isplane==0
        wc3=back_proj_line(centerx,kx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
        r3dx(1)= wc3(2);
        r3dx(2)= wc3(1);
        r3dx(3)= -wc3(3)+stdx;
    else
        wc3=back_proj_plane(centerx,kx,r3d0,stdx,sid,gAng,p2d(1),p2d(2));
        r3dx(1)= wc3(2);
        r3dx(2)= wc3(1);
        r3dx(3)= -wc3(3)+stdx;
    end    
    else % if the data in the database is not continuous, set the flag to 2.
    is_discon=2;
    r3dx(1)= 0;
    r3dx(2)= 0;
    r3dx(3)= 0;
end % of if (o3dx(jx+1,4)-o3d(jx,4)<0.2)
xx=range(ko3dd(1:kki-1,1:3),1);
range_o3d=xx(1,:);
return