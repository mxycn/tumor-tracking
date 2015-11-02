function temp_pt1=kv_imaging_new(jx,centerx,kk,y_min,y_max,imaging_angle,seg_err,sid,stdx,so3d)
% use kv imager to estimate the 3D position of jx-th point.
p2d=zeros(2,1);
wc3=[0 0 0 0];
temp_pt1=[0 0 0];
p2d(1:2)=proj([so3d(jx,1) so3d(jx,2) so3d(jx,3)],imaging_angle+90,sid,stdx);
p2d(1)=p2d(1) + seg_err * randn();
p2d(2)=p2d(2) + seg_err * randn();   
wc3=back_proj_line(centerx,kk,y_min,y_max,stdx,sid,imaging_angle+90,p2d(1),p2d(2));
temp_pt1(1)= wc3(2);
temp_pt1(2)= wc3(1);
temp_pt1(3)= -wc3(3)+stdx;
return
    
    
