% function [std_dev2,aberror2,ab_dev2,perct95_3d,perct95_2d,perct99_3d,perct99_2d] = simulation_segment(so3d,conditions)
function [std_dev2,aberror2,ab_dev2,parameters,threed,twod,WaperR_adjusted] = simulation_segment(so3d,conditions)
% simulate on one segment
% projecting and reconstructing the jth point.
%% initializing the parameters and arrays...
% is_pre = 1 if prediction is applied, 0 otherwise.
is_pre = conditions(1);

% is_block = 1 if block is simulated.
is_block = conditions(2);

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.
factor1 = conditions(3);

% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.
factor2 = conditions(4);

% factor3 specifies the ratio of the threshold controlling the predicted 
% changes of positions to the extent in their respective direction.
factor3 = conditions(5);

seg_err = conditions(6); % specifies the rms error of detection(=0.5mm).
lag_time = conditions(7); % the latent time
duration_ = conditions(10); % the time interval for one cycle of arc therapy (=72)
init_gAng = conditions(11); % initial angle(=-179)
end_gAng = conditions(12); % final angle(=+179)

% runlength: the number of initial points whose actual 3D coordinates
%  were measured in the first few seconds(~ 26)
init_time = conditions(13);
incrt = conditions(14);
runlength = floor(init_time/(incrt*0.0385))+1;
% one sample point are taken ever incrt data points, whose 
% period is approximately 0.038545s.incrt(= 4)
incrt = conditions(14);

sid = conditions(15); % sid: the Source-Imager distance.(=1500)
stdx = conditions(16); % stdx: the source-isocenter distance.(=1000)
r2threshold = conditions(23);
num_train = conditions(24);
marginx = conditions(25);

angular_speed=abs(end_gAng-init_gAng)/duration_; % the rate of gantry rotation 
ang_incrt=angular_speed*incrt*0.038545; % the angle interval between data points.
n_lag=round(lag_time/(incrt*38.5)); % the number of latent points
points_period=floor(358.0/ang_incrt);

count95_3 = 0;
count95_4 = 0;

r3d = zeros(points_period, 4);
r3dr = zeros(points_period, 4);
r3dp = zeros(points_period, 4);
incr_3d = zeros(points_period, 3);

kr = zeros(points_period,3); % the reconstructed direction vector 
ko = zeros(points_period,3); % the original direction vector
ym = zeros(points_period,2); % the end points of the PCA axis
centerx = zeros(points_period,3); % the center of the trajectory.
linearity = zeros(points_period,1);
is_kv = zeros(points_period,1);
% rerror3d = zeros(points_period,1);
% perror3d= zeros(points_period,1);
rerror3d = zeros(points_period,3);
perror3d= zeros(points_period,3);
perror2d= zeros(points_period,1);
is_plane = zeros(points_period,1);
predictionquality=zeros(points_period,1);
Bc = zeros(points_period,3);
cost = zeros(points_period,1);

aberror2=0;
ab_dev2=0;

T1=period_time(so3d(:,1),ang_incrt)
T2=period_time(so3d(:,2),ang_incrt)
T3=period_time(so3d(:,3),ang_incrt)
WaperR_adjusted=(periodic_evaluate(so3d(:,1),T1,1/ang_incrt)+...
                 periodic_evaluate(so3d(:,2),T2,1/ang_incrt)+...
                 periodic_evaluate(so3d(:,3),T3,1/ang_incrt))/3
%% Initializing
% add detection error to the points of the first few seconds, and calculate the 
% center point for the first few seconds.

%initializing the avarage
av3d=zeros(points_period,3);
av3d(1,:) = r3d(1,1:3);
if n_lag > 0.1
    r3dr(1:n_lag,:)=NaN;
    r3dp(1:n_lag,:)=NaN;
end
centerx(runlength,:)= mean(r3d(1:3,1));    
jj=1;
%% start simulate the position estimation and prediction point by point
while jj <= points_period
    % project who
    % start timing
    gAng=init_gAng+ang_incrt*(jj-1);    
    if jj<=runlength
        r3d(jj,1:3)=so3d(jj,1:3) + seg_err * stdx/sid * randn(1,3);
        
        if jj>2
            av3d(jj-1,:)=(r3d(jj-2,1:3)+2*r3d(jj-1,1:3)+r3d(jj,1:3))/4;  
        else
            av3d(1,:) = r3d(1,1:3);            
        end
        if jj > n_lag
            r3dr(jj,:)=r3d(jj-n_lag,:); 
            [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr); 
        end % predicting for the first few seconds
    else % of jj<=runlength. estimating position
        [coeff,score,roots] = princomp(r3d(1:runlength,1:3));
       dir_vector = coeff(:,1)';
       centerx0=mean(r3d(1:runlength,1:3),1);
       xfit0=repmat(centerx0(1,:),runlength,1) + score(:,1)*coeff(:,1)';
       for i=2:runlength
       if (xfit0(i,1)-xfit0(i-1,1))>=0
           r3d(i,4)=1;
       end
       end
%           so3d(jj,:)=so3d(jj-runlength,:);
        r3dx = [0 0 0 0];
        tic();
        [r3dx,centerx(jj-1,:),kr(jj-1,:),ko(jj-1,:),ym(jj-1,:),linearity(jj-1),is_plane(jj-1), Bc(jj,:),cost(jj)]= ...
            proj_backproj_new(seg_err,marginx,so3d,r3d,jj,gAng,runlength,sid,stdx);        
            xx1=toc();
%             timex(1)=timex(1)+xx1;
            if xx1>0.010
                fprintf('extra long estimation time:%f \n', xx1);
            end
        % accepting the returned values.
        r3d(jj,1)=r3dx(1);
        r3d(jj,2)=r3dx(2);
        r3d(jj,3)=r3dx(3);
        r3d(jj,4)=r3dx(4);
        
        % predicting the position the calculation is done starting from the j0+n_lag, the
        % first n_lag points are neglected.
        
       % smoothing the position
       av3d(jj-n_lag-1,:)=(r3d(jj-n_lag-2,1:3)+2*r3d(jj-n_lag-1,1:3)+r3d(jj-n_lag,1:3))/4;     
%         av3d(jj-n_lag,:)=(2*r3d(jj-n_lag-1,:)+3*r3d(jj-n_lag,:))/5;
            
      %% start to predict ..
        if n_lag>0
            tic();
            % calculate the immediate past range
            [extentx extentv]=near_range3d(jj,runlength,n_lag,r3d);                     
            % calculate the difference between the point n_log
            % before and the current point, for prediction. Note the
            % input information is from the points at least n_lag
            % before.
            [incr_3d(jj-n_lag,:) predictionquality(jj)]=predicting(is_pre,factor1,factor2,r2threshold,num_train,n_lag,jj-n_lag,extentx,extentv,av3d);    
            
            % Add the difference to the position n_lag points before 
            % to get the predicted position of the current point.
            % If the calculated difference is too large, replace
            % it with a threshold value, which is based on the
            % amplitude of motion in the first 4 seconds. 
            [r3dp(jj,:) predq]=prediction_screening(jj,n_lag,r3d,incr_3d,extentv,factor3);
            predictionquality(jj) = min(predictionquality(jj),predq);
            % optimize the prediction result with the updated 
            % correlation information n_lag points before 
            if jj>1+n_lag+runlength
                r3dr(jj,:)=prediction_optimizing_new(jj,r3dp,stdx,centerx(jj-n_lag,:),kr(jj-n_lag,:),linearity(jj-n_lag)); % true realtime  
            else
                r3dr(jj,:)=r3dp(jj,:);
            end
            
            % the following is the simulation of blockage and the
            % update of estimated coordinates and centerpoints.
            
            if is_block==1 
                if is_kv(jj+2-n_lag,1)==1 % that means kV image is available at jj+1-n_lag, 
                    % judge whether point jj-1 is actually blocked.                                           
                    if (jj+1-n_lag>202 && jj+1-n_lag<281) || (jj+1-n_lag>351 && jj+1-n_lag<430) 
                        temp_pt1=[0 0 0];
                        is_kv(jj+3-n_lag,1)=1; % more kV images to come
                        temp_pt1=kv_imaging_new(jj+1-n_lag,centerx(jj-n_lag,:),kr(jj-n_lag,:),ym(jj-n_lag,1),ym(jj-n_lag,2),gAng-ang_incrt,seg_err,sid,stdx,so3d);
                        r3d(jj+1-n_lag,:)=temp_pt1(:);
                    else % if point jj+1-n_lag is not blocked, stereostatic imaging is used.
                        r3d(jj+1-n_lag,1)=so3d(jj+1-n_lag,1) + seg_err * stdx/sid * randn();
                        r3d(jj+1-n_lag,2)=so3d(jj+1-n_lag,2) + seg_err * stdx/sid * randn();
                        r3d(jj+1-n_lag,3)=so3d(jj+1-n_lag,3) + seg_err * stdx/sid * randn();  
                        r3d(jj+1-n_lag,4)=so3d(jj+1-n_lag,4); 
                        is_kv(jj+2-n_lag,1)=0;                                 
                    end
                    % the point jj+1-n_lag needs to be updated.                
                    % see if jj+1-n_lag is blocked, if yes, calculate                
                    % its value, instead of measuring it.
                elseif (jj+1-n_lag>202 && jj+1-n_lag<281) || (jj+1-n_lag>351 && jj+1-n_lag<430) % of is_kv(jj+2-n_lag,1)==1 
                    r3d(jj+1-n_lag,1)=(r3dr(jj+1-n_lag,1)+2*r3d(jj-n_lag,1))/3.0;
                    r3d(jj+1-n_lag,2)=(r3dr(jj+1-n_lag,2)+2*r3d(jj-n_lag,2))/3.0;
                    r3d(jj+1-n_lag,3)=(r3dr(jj+1-n_lag,3)+2*r3d(jj-n_lag,3))/3.0;
                    is_kv(jj+4-n_lag,1)=1; % set current flag to 1, 
                    % meaning turn the kv imager on
                    % at jj(found at jj+1-n_lag, on at
                    % jj+2-n_lag, availabe at jj+3-n_lag).
                end % of is_kv(jj+2-n_lag,1)==1 
                % calculate the error.
%                 [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj-1,gAng-ang_incrt,sid,stdx,so3d,r3d,r3dr);  
                [rerror3d(jj,:) perror3d(jj,:) perror2d(jj)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
            else % of if is_block==1. No blockage simulation.  
%                 [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
                [rerror3d(jj,:) perror3d(jj,:) perror2d(jj)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
            end % of if is_block==1. blockage simulation is finished. 
            xx2=toc();
%             timex(2)=timex(2)+xx2;
            if xx2>0.02
                fprintf('extra long prediction time:%f \n', xx2);
            end
        else % of n_lag>0. No latency is considered
            % if the point position is measured in the first 4
            % seconds, or no latency is considered.
            % no prediction is made, use the latest point estimated  
            r3dr(jj,1:3)=r3d(jj,1:3);
            % calculate the projection of estimated position.
%             [rerror3d(jj) perror3d(jj) perror2d(jj)]=error_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
            [rerror3d(jj,:) perror3d(jj,:) perror2d(jj)]=errorv_calculation(jj,gAng,sid,stdx,so3d,r3d,r3dr);
        end % of n_lag>0.
    end % of jj<=runlength
    %         fprintf('extra long prediction time22:%f %u  \n',  r3dr(jj,1), jj);
    aberror2=aberror2+perror3d(jj);
    ab_dev2=ab_dev2+perror2d(jj);   
    if jj == 1 + runlength
        threed = rerror3d(jj,:);
        twod = perror3d(jj,:);      
%         threed = sqrt(rerror3d(jj));
%         twod = sqrt(perror2d(jj));   
        speed = norm(r3dr(jj-1,:)-r3dr(jj-3,:));
        cos1 = Bc(jj,:)*kr(jj-2,:)'; % cosine of the angle between the principal axis and the treatment beam
%         cos2 = Bc(jj,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed; % cosine of the angle between the velocity axis and the beam
        cos1s = cos1*cos1; 
%         cos2s = cos2*cos2;
        preerr = sqrt(perror2d(jj-2));        
        parameters = [speed cos1s preerr predictionquality(jj),cost(jj-2)];
    elseif jj > 1 + runlength
%          elseif jj > 1 + runlength && mod(jj,2) == 0
        threed = [threed;rerror3d(jj,:)];
        twod = [twod;perror3d(jj,:)];
        speed = norm(r3dr(jj-1,:)-r3dr(jj-3,:));
        cos1 = Bc(jj,:)*kr(jj-2,:)';
%         cos2 = kr(jj-1,:)*(r3dr(jj,:)-r3dr(jj-1,:))'/speed;
        cos1s = cos1*cos1;
%         cos2s = cos2*cos2;
        preerr = sqrt(perror2d(jj-2));
        parameters = [parameters;[speed cos1s preerr predictionquality(jj),cost(jj-2)]];        
    end
    jj=jj+1;
    %     denomatorx=sqrt(ko(jj-2,1)*ko(jj-2,1)+ko(jj-2,2)*ko(jj-2,2)+ko(jj-2,3)*ko(jj-2,3))*...
    %         sqrt(ko1*ko1+ko2*ko2+ko3*ko3);
    %     cosangle=(ko(jj-2,1)*ko1+ko(jj-2,2)*ko2+ko(jj-2,3)*ko3)/denomatorx;
    %     if abs(cosangle)<0.966
    %         num_coord=num_coord+1;                  
    %     end
end % of "while jj< = point_period."
% calculate the variation of the trajectory itself.
od=zeros(points_period,3);
% construct a new array for calculating the standard deviation
% of the trajectory.
od(:,:)=so3d(1:points_period,:);    
% sum up the square of all standard deviation.
variancex=std(od(1:points_period,1))*std(od(1:points_period,1))+...
    std(od(1:points_period,2))*std(od(1:points_period,2))+...
    std(od(1:points_period,3))*std(od(1:points_period,3));
std_dev2=variancex; 
aberror2 = aberror2/(points_period-n_lag);
ab_dev2 = ab_dev2/(points_period-n_lag);

%% plotting
xcoord=zeros(points_period,1);
for i=1:points_period
    xcoord(i)=ang_incrt*(i-1)+init_gAng+180;
end

% std1=norm(perror2d(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); % average 2D error with prediction
% max1=max(perror2d(j0+runlength+2:j0+points_period-1));
% std2=norm(perror3dx(j0+n_lag:j0+points_period-1))/sqrt(points_period-n_lag); 
% std3=norm(recon_error3d(j0+runlength:j0+points_period-1))/sqrt(points_period-runlength); % average 3D error of pure estimation

perct95_3d=sqrt(get_percentile_error(perror3d,0.95,points_period));
perct95_2d=sqrt(get_percentile_error(perror2d,0.95,points_period));
perct99_3d=sqrt(get_percentile_error(perror3d,0.99,points_period));
perct99_2d=sqrt(get_percentile_error(perror2d,0.99,points_period));

% if percent95_2dx>3.0
%     count95_3=count95_3+1;
%     if percent95_2dx>4.0
%         count95_4=count95_4+1;
%     end
% end
%              plotting_condition=((std1 > 1.35 & std2 < 2.45) | (std1 > 12.69 & std2>13.29) | max1 > 19.3) &...
%                  ((kk1(jj-3,4)-kk1(jj-3,3))>3.0*(kkk1(4)-kkk1(3)) | ...
%                  (kk1(jj-3,8)-kk1(jj-3,7))>3.0*(kkk1(8)-kkk1(7)) | ...
%                  (kk2(jj-3,8)-kk2(jj-3,7))>3.0*(kkk2(8)-kkk2(7))) & ...
%                  (abs(kko1(jj-1,9))>0.92 | abs(kko2(jj-1,9))< 0.92 | ...
%                  (abs(kko1(jj-1,9))>0.9 & abs(kko2(jj-1,9))> 0.9));

%              conditionx=std1<1.70 & std_deviation> 3.0 & k==1 & (j0==1871 | j0==11774);  
conditionx= variancex > 0.30;  
% %              conditionx= k==1 & j0<500;
i0= runlength+1;
modex = 1;
n_segment = 1;
% plotfigures(conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,sqrt(perror2d),sqrt(perror3d),sqrt(rerror3d));
plotfigures(i0,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,[zeros(i0,1);parameters(:,2)],sqrt(perror3d),sqrt(rerror3d),runlength,ang_incrt,centerx,kr);
modex = 2;
n_segment = 2;
% plotfigures(conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,sqrt(perror2d),sqrt(perror3d),sqrt(rerror3d));
plotfigures(i0,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,[zeros(i0,1);parameters(:,2)],sqrt(perror3d),sqrt(rerror3d),runlength,ang_incrt,centerx,kr);
pause;


    
        
    