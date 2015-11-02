function [extentx extentv]=near_range3d(jx,runlength,n_lag,r3d)
% calculate the immediate past range of the points.
% Input: jx-the current point; runlength-the number of points in the first 4 seconds;
%        n_lag: latency points; r3d: estimated 3D coordinates
% if jx>j0+16*runlength
if jx>1+2*runlength
%     starti=jx-16*runlength;
    starti=jx-2*runlength;
else
    starti=2;
end
% calculate the maximum velocity along the 3 directions.
min3d=[min(r3d(starti:jx-n_lag-1,1)) min(r3d(starti:jx-n_lag-1,2)) min(r3d(starti:jx-n_lag-1,3))];
max3d=[max(r3d(starti:jx-n_lag-1,1)) max(r3d(starti:jx-n_lag-1,2)) max(r3d(starti:jx-n_lag-1,3))];
% min3d=[min(r3d(starti:starti+runlength-1,1)) min(r3d(starti:starti+runlength-1,2)) min(r3d(starti:starti+runlength-1,3))];
% max3d=[max(r3d(starti:starti+runlength-1,1)) max(r3d(starti:starti+runlength-1,2)) max(r3d(starti:starti+runlength-1,3))];
extentx=[max3d(1)-min3d(1) max3d(2)-min3d(2) max3d(3)-min3d(3)]; 

extentv = max(abs(r3d(starti:jx-n_lag-1,:) - r3d(starti-1:jx-n_lag-2,:)),[],1);
%     max(abs(r3d(starti-1:jx-n_lag-1,2) - r3d(starti-1:jx-n_lag-2,2))) ...
%     max(abs(r3d(starti-1:jx-n_lag-1,3) - r3d(starti-1:jx-n_lag-2,3)))];
% pause(1);