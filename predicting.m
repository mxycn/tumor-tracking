function [incrx predictionquality]=predicting(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jx,extentx,extentv,av3d)
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

% the regression equation is
%     incrx(:)=a+b*av3d(jx-1,:)+c*(av3d(jx-3,:);
        
incrx = [0 0 0];
predictionquality = 1;
if is_pre == 1
    delta_y_theta_k =[0 0 0];
    delta_y_theta_k_l = [0 0 0];
    delta_y_theta_k_2l = [0 0 0];
    delta_y_theta_k(:)= av3d(jx-n_lag,:)-av3d(jx,:);
    delta_y_theta_k_l(:)=av3d(jx-n_lag-1,:)-av3d(jx-1,:);
    delta_y_theta_k_2l(:)=av3d(jx-1-n_lag-1,:)-av3d(jx-2,:);
    
    % The variables entering the regression, selected from the raw data,
    % the last being the dependent variable.
    v3d1=zeros(jx,3);
    v3d2=zeros(jx,3);
    v3d=zeros(jx,3);
    % The indices for different points in the regression.
    k=[0 0 0];    
    % Construct the independent variables and the dependent variable for the 
    % regression in the SI direction.
    % Consider only those points that are similar to the present point in terms
    % of position and velocity.
    % ii is the index for the input data, which is screened by the conditions.
    % Because the regression in a sense is an autoregression, and the input data
    % involves av3d(ii-3,1), and av3d(ii+n_lag,1), ii runs from 4 to jx-1-n_lag-3,
    % where -3 is added to eliminate the points close to the current point.
    
    % See file "details" for explanation, here av3d(ii-1,1)+av3d(ii-2,1) replaces
    % av3d(ii-1,1).
    if extentx(1) < 20       
    for i3 = 1:3  
        mark_vm=[0 0 0];
        for ii=4:jx-1-n_lag
            % the reference for the latest velocity        
            mark_vm(i3)=av3d(jx-1,i3)-av3d(jx-2,i3);  
            % av3d(ii-1,1) plays the same role in the regression as what av3d(jx-1,1)
            % plays in the output. The first condition compares the position, the 
            % second compares the velocity
            if abs(av3d(ii-1,i3)-av3d(jx-1,i3))<extentx(i3)/factor1 && ...
                    abs(av3d(ii-1,i3)-av3d(ii-2,i3)-mark_vm(i3))<extentv(i3)*factor2 
                k(i3)=k(i3)+1;
                v3d1(k(i3),i3)=av3d(ii-1,i3); % the first independent variable.
                v3d2(k(i3),i3)=av3d(ii-3,i3); % the second independent variable.
                v3d(k(i3),i3)=av3d(ii+n_lag,i3)-av3d(ii,i3);   % the dependent variable.
            end
        end    
        % two-variate regression
        if k(i3)>num_train % if more than 3 set of numbers are acquired as the input data, do regression            
            [b r2]=two_var_regress(k(i3),v3d1(1:k(i3),i3),v3d2(1:k(i3),i3),v3d(1:k(i3),i3));  
            if r2<predictionquality
                predictionquality = r2;
            end
            if r2>r2threshold
                incrx(i3)=b(1)+b(2)*av3d(jx-1,i3)+b(3)*av3d(jx-3,i3);                   
            else
                 incrx(i3)=mean(v3d(1:k(i3),i3));                  
            end            
        elseif k(i3)>0.2 % if only one, two or three sets of numbers are acuquired, use the average
            incrx(i3)=mean(v3d(1:k(i3),i3));  
            predictionquality = 0;
%             if no reference point, predict according to the shape(very rough, MULIN algorithm)
        else
            incrx(i3)=0;
            predictionquality = 0;
        end
    end
    else % use only the motion in SI direction as the independent vairable
        k1 = 0;
        for ii=4:jx-1-n_lag
            % the reference for the latest velocity        
            mark_vm(1)=av3d(jx-1,1)-av3d(jx-2,1);  
            % av3d(ii-1,1) plays the same role in the regression as what av3d(jx-1,1)
            % plays in the output. The first condition compares the position, the 
            % second compares the velocity
            if abs(av3d(ii-1,1)-av3d(jx-1,1))<extentx(1)/factor1 && ...
                    abs(av3d(ii-1,1)-av3d(ii-2,1)-mark_vm(1))<extentv(1)*factor2 
                k1=k1+1;
                v3d1(k1,1)=av3d(ii-1,1); % the first independent variable.
                v3d2(k1,1)=av3d(ii-3,1); % the second independent variable.
                v3d(k1,:)=av3d(ii+n_lag,:)-av3d(ii,:);   % the dependent variable.                
            end
        end    
        % two-variate regression
        if k1>num_train % if more than 3 set of numbers are acquired as the input data, do regression
            [b1 r21]=two_var_regress(k1,v3d1(1:k1,1),v3d2(1:k1,1),v3d(1:k1,1));  
            [b2 r22]=two_var_regress(k1,v3d1(1:k1,1),v3d2(1:k1,1),v3d(1:k1,2));  
            [b3 r23]=two_var_regress(k1,v3d1(1:k1,1),v3d2(1:k1,1),v3d(1:k1,3));  
            predictionquality = (r21+r22+r23)/3;            
            if r21>r2threshold
                incrx(1)=b1(1)+b1(2)*av3d(jx-1,1)+b1(3)*av3d(jx-3,1); 
            else
                 incrx(1)=mean(v3d(1:k1,1)); 
            end            
            if r22>r2threshold
                incrx(2)=b2(1)+b2(2)*av3d(jx-1,1)+b2(3)*av3d(jx-3,1); 
            else
                 incrx(2)=mean(v3d(1:k1,2)); 
            end            
            if r23>r2threshold
                incrx(3)=b3(1)+b3(2)*av3d(jx-1,1)+b3(3)*av3d(jx-3,1); 
            else
                 incrx(3)=mean(v3d(1:k1,3)); 
            end            
        elseif k1>0.2 % if only one, two or three sets of numbers are acuquired, use the average
            incrx(:)=mean(v3d(1:k1,:));    
            predictionquality = 0;
        else
            incrx(:)=0;
            predictionquality = 0;
        end
    end
end