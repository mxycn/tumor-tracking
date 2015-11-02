function [incr1 incr2 incr3 isfail]=predicting(is_pre,factor1,factor2,n_lag,jx,j0,extentx,av3d)
% This function returns the predicted changes in the three directions of 
% motion, incr1, incr2, and incr3, which will be added to the coordinates 
% n_lag points earlier in function workflow(...).

% is_pre = 1 if prediction is applied, 0 otherwise.

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.

% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.

% n_lag is the number of latent points, n=3 corresponds to 460ms, n=2
% cooresponds to 310ms.

% extentx gives the range of motion in the first 4 seconds.

% The following are coefficients of regression equation, for the SI, LR, AP
% directions. Two-variate regression is better than three-variate
% regression, so only three coefficient for each direction of motion is
% neccessary.
a1=0;
b1=0;
c1=0;
d1=0;
e1=0;
f1=0;
a2=0;
b2=0;
c2=0;
d2=0;
e2=0;
f2=0;
a3=0;
b3=0;
c3=0;
d3=0;
e3=0;
f3=0;
isfail=1;
delta_y_theta_k=zeros(3,1);
delta_y_theta_k_l=zeros(3,1);
delta_y_theta_k_2l=zeros(3,1);
delta_y_theta_k(1)=av3d(jx-n_lag-j0+1,1)-av3d(jx-j0+1,1);
delta_y_theta_k(2)=av3d(jx-n_lag-j0+1,2)-av3d(jx-j0+1,2);
delta_y_theta_k(3)=av3d(jx-n_lag-j0+1,3)-av3d(jx-j0+1,3);
delta_y_theta_k_l(1)=av3d(jx-n_lag-j0,1)-av3d(jx-j0,1);
delta_y_theta_k_l(2)=av3d(jx-n_lag-j0,2)-av3d(jx-j0,2);
delta_y_theta_k_l(3)=av3d(jx-n_lag-j0,3)-av3d(jx-j0,3);
delta_y_theta_k_2l(1)=av3d(jx-1-n_lag-j0,1)-av3d(jx-1-j0,1);
delta_y_theta_k_2l(2)=av3d(jx-1-n_lag-j0,2)-av3d(jx-1-j0,2);
delta_y_theta_k_2l(3)=av3d(jx-1-n_lag-j0,3)-av3d(jx-1-j0,3);

% The averaged coordinates, initialization. 

% The variables entering the regression, the last being the dependent
% variable, initialization.
v3d1=zeros(jx-j0+1,3);
v3d2=zeros(jx-j0+1,3);
v3d=zeros(jx-j0+1,3);
% The indices for different points in the regression.
k1=0;
k2=0;
k3=0;
% Initializing the returns.
incr1=0;
incr2=0;
incr3=0;

% Construct the independent variables and the dependent variable for the 
% regression in the SI direction.
% Consider only those points that are similar to the present point in terms
% of position and velocity.
% ii is the index for the input data, which is screened by the conditions.
% Because the regression in a sense is an autoregression, and the input data
% involves av3d(ii-3,1), and av3d(ii+n_lag,1), ii runs from 4 to jx-j0-n_lag-3,
% where -3 is added to eliminate the points close to the current point.

% See file "details" for explanation, here av3d(ii-1,1)+av3d(ii-2,1) replaces
% av3d(ii-1,1).
mark_sm1=av3d(jx-j0,1)-av3d(jx-j0-1,1);     % a reference for the latest averaged velocity
for ii=4:jx-j0-n_lag
    % av3d(ii-1,1) plays the same role in the regression as what av3d(jx-j0,1)
    % plays in the output. The first condition compares the position, the 
    % second compares the velocity
    if abs(av3d(ii-1,1)-av3d(jx-j0,1))<extentx(1)/factor1 & ...
            abs(av3d(ii-1,1)-av3d(ii-2,1)-mark_sm1)<extentx(1)/factor2 & ...
        ~isnan(av3d(ii+n_lag,1)-av3d(ii,1)) & ~isnan(av3d(ii-3,1))
        k1=k1+1;
        v3d1(k1,1)=av3d(ii-1,1); % the first independent variable.
        v3d2(k1,1)=av3d(ii-1,1)-av3d(ii-3,1); % the second independent variable.
        v3d(k1,1)=av3d(ii+n_lag,1)-av3d(ii,1);   % the dependent variable.
    end
end

% two-variate regression
if k1>4.1 & ~isnan(av3d(jx-j0-2,1)) % if more than 3 set of numbers are acquired as the input data, do regression
    [a1 b1 c1]=two_var_regression(k1,v3d1(1:k1,1),v3d2(1:k1,1),v3d(1:k1,1));  
    incr1=a1+b1*av3d(jx-j0,1)+c1*(av3d(jx-j0,1)-av3d(jx-j0-2,1));
elseif k1>0.2 % if only one, two or three sets of numbers are acuquired, use the average
    incr1=mean(v3d(1:k1,1));    
% if no reference point, predict according to the shape(very rough, MULIN algorithm)
else
    delta_y_theta_k(1)=av3d(jx-n_lag-j0+1,1)-av3d(jx-j0+1,1);
    delta_y_theta_k_l(1)=av3d(jx-n_lag-j0,1)-av3d(jx-j0,1);
    delta_y_theta_k_2l(1)=av3d(jx-1-n_lag-j0,1)-av3d(jx-1-j0,1);
    incr1=-2*delta_y_theta_k(1)+delta_y_theta_k_l(1);
    isfail=2;
    if isnan(incr1)
        incr1=0;
    end
end

% Construct the independent variables and the dependent variable for the 
% regression in the LR direction.
mark_sm1=av3d(jx-j0,2)-av3d(jx-j0-1,2);
for ii=4:jx-j0-n_lag
    if abs(av3d(ii-1,2)-av3d(jx-j0,2))<extentx(2)/factor1 & ...
            abs(av3d(ii-1,2)-av3d(ii-2,2)-mark_sm1)<extentx(2)/factor2 & ...
            ~isnan(av3d(ii+n_lag,2)-av3d(ii,2)) & ~isnan(av3d(ii-3,2))
        k2=k2+1;
        v3d1(k2,2)=av3d(ii-1,2);
        v3d2(k2,2)=av3d(ii-1,2)-av3d(ii-3,2);
        v3d(k2,2)=av3d(ii+n_lag,2)-av3d(ii,2);
    end
end
 
% two-variate regression
if k2>4.1 & ~isnan(av3d(jx-j0-2,2))
    [a2 b2 c2]=two_var_regression(k2,v3d1(1:k2,2),v3d2(1:k2,2),v3d(1:k2,2));
    incr2=a2+b2*av3d(jx-j0,2)+c2*(av3d(jx-j0,2)-av3d(jx-j0-2,2));
elseif k2>0.2
    incr2=mean(v3d(1:k2,2));
else
    delta_y_theta_k(2)=av3d(jx-n_lag-j0+1,2)-av3d(jx-j0+1,2);
    delta_y_theta_k_l(2)=av3d(jx-n_lag-j0,2)-av3d(jx-j0,2);
    delta_y_theta_k_2l(2)=av3d(jx-1-n_lag-j0,2)-av3d(jx-1-j0,2);
    incr2=-2*delta_y_theta_k(2)+delta_y_theta_k_l(2);
    isfail=2;
    if isnan(incr2)
        incr2=0;
    end
end

% Construct the independent variables and the dependent variable for the 
% regression in the AP direction.
mark_sm1=av3d(jx-j0,3)-av3d(jx-j0-1,3);
for ii=4:jx-j0-n_lag
    if abs(av3d(ii-1,3)-av3d(jx-j0,3))<extentx(3)/factor1 & ...
            abs(av3d(ii-1,3)-av3d(ii-2,3)-mark_sm1)<extentx(3)/factor2 & ...
            ~isnan(av3d(ii+n_lag,3)-av3d(ii,3)) & ~isnan(av3d(ii-3,3))
        k3=k3+1;
        v3d1(k3,3)=av3d(ii-1,3);
        v3d2(k3,3)=av3d(ii-1,3)-av3d(ii-3,3);
        v3d(k3,3)=av3d(ii+n_lag,3)-av3d(ii,3);
    end
end

% two-variate regression
if k3>4.1 & ~isnan(av3d(jx-j0-2,3))
    tic();
    [a3 b3 c3]=two_var_regression(k3,v3d1(1:k3,3),v3d2(1:k3,3),v3d(1:k3,3));
    xxx=toc();
    if xxx>0.006
        fprintf('large regression time observed: %f \n',xxx);
    end        
    incr3=a3+b3*av3d(jx-j0,3)+c3*(av3d(jx-j0,3)-av3d(jx-j0-2,3));
elseif k3>0.2
    incr3=mean(v3d(1:k3,3));
else
    delta_y_theta_k(3)=av3d(jx-n_lag-j0+1,3)-av3d(jx-j0+1,3);
    delta_y_theta_k_l(3)=av3d(jx-n_lag-j0,3)-av3d(jx-j0,3);
    delta_y_theta_k_2l(3)=av3d(jx-1-n_lag-j0,3)-av3d(jx-1-j0,3);
    incr3=-2*delta_y_theta_k(3)+delta_y_theta_k_l(3);
    isfail=2;
    if isnan(incr3)
        incr3=0;
    end
end

% If no prediction is made, the changes are set to zero.
if is_pre<0.5
    incr1=0;
    incr2=0;
    incr3=0;
end
return