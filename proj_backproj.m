function [r3dx,is_discon,kko1,kko2,kk1,kk2,is_jump,iskv]=proj_backproj(seg_err,...
    k,o3d,r3d,jx,j0,init_gAng,end_gAng,delta_angle,ang_incrt,runlength,sid,stdx);
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
is_jump=0;
iskv=0;
% the current gantry angle.
gAng=init_gAng+ang_incrt*(jx-j0);
r3d0=zeros(jx-j0+1,3);

% the following stores real-time correlation information in the original 
% 3D data, used to screen trajectories.
kko1=[0 0 0 0 0 0 0 0 0];
kko2=[0 0 0 0 0 0 0 0 0];
kko1=find_corr_info(jx-j0,o3d(k,j0:jx-1,1),...
    o3d(k,j0:jx-1,2)); % in left + (data coordinate y, liuwu coordinate x)
kko2=find_corr_info(jx-j0,o3d(k,j0:jx-1,1),...
    stdx-o3d(k,j0:jx-1,3));

% the following stores real-time correlation information in the estimated 
% 3D data, used to update the correlation informaiton.
kk1=[0 0 0 0 0 0 0 0 0];
kk2=[0 0 0 0 0 0 0 0 0];
kk3=[0 0 0 0 0 0 0 0 0];

% if the velocity is too large, abandon the current trial
if abs(o3d(k,jx,1)-o3d(k,jx-1,1)) > 6.5 | ...
        abs(o3d(k,jx,2)-o3d(k,jx-1,2)) > 6.5 | ...
        abs(o3d(k,jx,3)-o3d(k,jx-1,3)) > 6.5
    is_jump=1; 
    is_discon=2;
    r3dx(1)= 0;
    r3dx(2)= 0;
    r3dx(3)= 0;
    return
end

for ix=1:runlength    
    r3d0(ix,1)=r3d(j0+ix-1,2);
    r3d0(ix,2)=r3d(j0+ix-1,1);
    r3d0(ix,3)=stdx-r3d(j0+ix-1,3);
end

p2d=[0 0];
wc3=[0 0 0 0];
% if the point data in the database is continuous(based on the time 
% information), perform estimation, flag it otherwise.
if (o3d(k,jx+1,4)-o3d(k,jx,4)<0.2) & (gAng<=end_gAng) 
    o_3d=[0 0 0 0];
    o_3d(:)=o3d(k,jx,:);
    p2d(1:2)=proj([o_3d(1) o_3d(2) o_3d(3)],gAng,sid,stdx);
    % kk can be updated with special considerations. delta_angle=20 when duration=72.
    kk1=find_corr_info(jx-j0,r3d(j0:jx-1,1),r3d(j0:jx-1,2)); % in left + (data coordinate y, liuwu coordinate x)
    kk2=find_corr_info(jx-j0,r3d(j0:jx-1,1),stdx-r3d(j0:jx-1,3)); % Anterior+: data coordinate z+, Liuwu coordinate z-        
%     
%     kk1=find_corr_info(runlength,r3d(j0:runlength+j0-1,1),r3d(j0:runlength+j0-1,2)); % in left + (data coordinate y, liuwu coordinate x)
%     kk2=find_corr_info(runlength,r3d(j0:runlength+j0-1,1),stdx-r3d(j0:runlength+j0-1,3)); % Anterior+: data coordinate z+, Liuwu coordinate z-    
    %% add segmentation error to the data
    p2d(1)=p2d(1) + seg_err * randn();
    p2d(2)=p2d(2) + seg_err * randn();
    wc3=back_proj(kk1,kk2,stdx,sid,gAng,p2d(1),p2d(2));
    if wc3(4)>0.9
        iskv=2;
        r3dx(1)=o3d(k,jx,1) + seg_err * stdx/sid * randn();
        r3dx(2)=o3d(k,jx,2) + seg_err * stdx/sid * randn();
        r3dx(3)=o3d(k,jx,3) + seg_err * stdx/sid * randn();
    else
        r3dx(1)= wc3(2);
        r3dx(2)= wc3(1);
        r3dx(3)= -wc3(3)+stdx;
    end        
%     if seg_err<0.41 & abs(kk1(9)*kk2(9))<0.49 % the following method works well only when seg_err is small.         
%         wc3=back_proj_plane(r3d0,stdx,sid,gAng,p2d(1),p2d(2));
%         r3dx(1)= wc3(2);
%         r3dx(2)= wc3(1);
%         r3dx(3)= -wc3(3)+stdx; 
%     else
%        wc3=back_proj(kk1,kk2,stdx,sid,gAng,p2d(1),p2d(2));
%         r3dx(1)= wc3(2);
%         r3dx(2)= wc3(1);
%         r3dx(3)= -wc3(3)+stdx;
%     end    
else % if the data in the database is not continuous, set the flag to 2.
    is_discon=2;
    r3dx(1)= 0;
    r3dx(2)= 0;
    r3dx(3)= 0;
    kk1=[0 0 0 0 0 0 0 0 0];
    kk2=[0 0 0 0 0 0 0 0 0];    
end % of if (o3dx(jx+1,4)-o3d(jx,4)<0.2)
return

