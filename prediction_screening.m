function [r3dp predictionquality]=prediction_screening(jx,n_lag,r3d,incr_3d,extentv,factor3)
% to screen out extraordinary predictions
% input: jx-the no. of current point; n_lag- the number of lag points; r3d:
%        the 3D coordinates on which the prediction is based; extentx:
%        reference range of variation; factor3: a factor to adjust the
%        ratio of the predicted change to the reference range. incr_3d: the
%        predicted change prior to the screening.
% output: the predicted 3D coordinates.

% SI motion
r3dp=[0 0 0 r3d(jx-n_lag,4)];
predictionquality = 1;
for i3 = 1:3
    maxi(i3)=extentv(i3)*factor3;
    if incr_3d(jx-n_lag,i3)>maxi(i3)
        r3dp(i3)=r3d(jx-n_lag,i3)+maxi(i3);
        if r3d(jx-n_lag,i3)-r3d(jx-2*n_lag)<-extentv(i3)*0.6 %If the previous velocity is small, then threshold acceleration is exceeded
            predictionquality = 0.3;
        end
    elseif incr_3d(jx-n_lag,i3)<-maxi(i3);
        r3dp(i3)=r3d(jx-n_lag,i3)-maxi(i3);
        if r3d(jx-n_lag,i3)-r3d(jx-2*n_lag)>extentv(i3)*0.6 %If the previous velocity is large, then threshold acceleration is exceeded
            predictionquality = 0.3;
        end
    else
        r3dp(i3)=r3d(jx-n_lag,i3)+incr_3d(jx-n_lag,i3);
        if abs(incr_3d(jx-n_lag,i3) + r3d(jx-n_lag,i3)-r3d(jx-2*n_lag)) > extentv(i3)*factor3*factor3           
            predictionquality = 0.3;
        end
    end
end