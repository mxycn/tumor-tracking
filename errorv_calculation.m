function [estimation_error prediction_error projected_error]=errorv_calculation(jx,g_Ang,sid,stdx,o3d,r3d,r3dr)
% calculate the estimation error(3D), the prediction error(2D), and the
% 2D beam's eye view error.
% Input: jx-the current no. of point; g_Ang-current gantry angle; 
%        sid-source-to-imager distance; o3d-the true 3D coordinates; 
%        r3d- the estimated 3D coordinates; r3dr- the predicted 3D
%        coordinates
% Output: estimation error, prediction error, projected error
pr_2d=[0 0 0];
po_2d=[0 0 0];
% calcualte the projected coordinates on the imaging plate.
po_2d(1:3)=proj([o3d(jx,1) o3d(jx,2) o3d(jx,3)],g_Ang,sid,stdx); 
pr_2d(1:3)=proj([r3dr(jx,1) r3dr(jx,2) r3dr(jx,3)],g_Ang,sid,stdx); 
% calculate the 3D error of the current point(reconstructed)
estimation_error = [r3d(jx,1)-o3d(jx,1) r3d(jx,2)-o3d(jx,2) r3d(jx,3)-o3d(jx,3)];
% calculate the 3D error of the current point(predicted)
prediction_error = [r3dr(jx,1)-o3d(jx,1) r3dr(jx,2)-o3d(jx,2) r3dr(jx,3)-o3d(jx,3)];
projected_error =((pr_2d(1)-po_2d(1))^2+...
    (pr_2d(2)-po_2d(2))^2)/(pr_2d(3)*pr_2d(3)); 

