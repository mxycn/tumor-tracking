function [r3dpx r3dpy r3dpz]=prediction_screening(jx,n_lag,r3d,incr_3d,extentx,factor3)
% to screen out extraordinary predictions
% input: jx-the no. of current point; n_lag- the number of lag points; r3d:
%        the 3D coordinates on which the prediction is based; extentx:
%        reference range of variation; factor3: a factor to adjust the
%        ratio of the predicted change to the reference range. incr_3d: the
%        predicted change prior to the screening.
% output: the predicted 3D coordinates.

% SI motion
maxi(1)=extentx(1)/factor3;
if incr_3d(jx-n_lag,1)>maxi(1)
    r3dpx=r3d(jx-n_lag,1)+maxi(1);
elseif incr_3d(jx-n_lag,1)<-maxi(1);
    r3dpx=r3d(jx-n_lag,1)-maxi(1);
else
    r3dpx=r3d(jx-n_lag,1)+incr_3d(jx-n_lag,1);
end

% LR motion
maxi(2)=extentx(2)/factor3;
if incr_3d(jx-n_lag,2)>maxi(2)
    r3dpy=r3d(jx-n_lag,2)+maxi(2);
elseif incr_3d(jx-n_lag,2)<-maxi(2)
    r3dpy=r3d(jx-n_lag,2)-maxi(2);
else
    r3dpy=r3d(jx-n_lag,2)+incr_3d(jx-n_lag,2);
end

% AP motion
maxi(3)=extentx(3)/factor3;
if incr_3d(jx-n_lag,3)>maxi(3)
    r3dpz=r3d(jx-n_lag,3)+maxi(3);
elseif incr_3d(jx-n_lag,3)<-maxi(3)
    r3dpz=r3d(jx-n_lag,3)-maxi(3);
else
    r3dpz=r3d(jx-n_lag,3)+incr_3d(jx-n_lag,3);
end