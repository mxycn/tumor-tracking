function [ndir,nactual,n_seg,rms,std_dev,devp,pct95_3d,pct95_2d,pct99_3d,pct99_2d,count95_3,timex] = ...
    simulation_fixed(is_pre,is_block,factor1,factor2,factor3,seg_err,lag_time,linearity,...
    amp,runlength,incrt,sid,stdx,angleset,onofftime,path_in)
% This function performs data input, reconstruction, and prediction, 
% error calculation and plotting

% is_pre = 1 if prediction is applied, 0 otherwise.

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.
% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.
% factor3 specifies the ratio of the threshold controlling the predicted 
% changes of positions to the extent in their respective direction.

% seg_err specifies the rms error of detection(=0.5mm).

% n_lag is the number of latent points, n=3 corresponds to 460ms, n=2
% cooresponds to 310ms.

% linearity presents the threshold of correlation coefficient for the motion of 
% the first 4 seconds. If the threshould is surpassed, the data of the 4
% second is discarded. The program moves to the data of next 4 seconds.

% amp specifies the maximum ratio of motion amplitude of later data points
% to that of the first 4 seconds.

% runlength(= 26): the number of initial points whose actual 3D coordinates
%                  were measured in the first 4 seconds;
% incrt(= 4): one sample point are taken ever incrt data points, whose 
%             period is approximately 0.038545s.
% sid(=1500): the Source-Imager distance.
% stdx(=1000): the source-isocenter distance.
% path_in: the folder of the data
n_lag=round(lag_time/(incrt*38.5));
% est_length=floor(3.0*runlength/5.0)-4;
 % the number of points used
nactual=0; % the number of points available
 % the sum of the square of standard deviation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModelerDataDir = path_in; % the folder of data
cd(ModelerDataDir); 
dbfolders = dir('DB*');
ndir = length(dbfolders); % the number of folders
npts=zeros(ndir,1); % the number of points in each folder
% define an array to accept the data
o3d=zeros(ndir,300000,4); % the array storing the original point data
npoints=zeros(ndir,1);
aberror2=zeros(ndir,1); % the sum of the square of 3D error
ab_dev2=zeros(ndir,1);  % the sum of the square of 2D error
std_dev2=zeros(ndir,1);
percent95_3d=zeros(ndir,1);
percent95_2d=zeros(ndir,1);
percent99_3d=zeros(ndir,1);
percent99_2d=zeros(ndir,1);
% read data from the data file
[o3d,npts]=read_data(ndir,incrt,path_in);
count95_3=0;
count95_4=0;
timex=[0 0];
num_plot=0;
% where npts is the number of points which equals the total original data points
% divided by incrt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_seg=0;                  % the number of trajectories
totaltime=onofftime(length(onofftime));
transition_time=50;
points_period=floor(totaltime/4/0.0385);
points_transition=round(transition_time/4/0.0385);
% k is the number of fractions
% for k=1:1  % ndir is the total number of fractions for that patient.
for k=1:ndir
    pr_2d=zeros(npts(k),2); % the projection of the assumed point after repositioning.
    r3d=zeros(npts(k),4);   % the 3D position estimated.
    pltr3d=zeros(npts(k),4);   % the 3D position estimated.
    plto3d=zeros(npts(k),4);   % the 3D position estimated.
    incr_3d=zeros(npts(k),3); % the predicted difference used for compensate for the lag of time.
    r3dp=zeros(npts(k),3); % the predicted 3D position without optimization.
    r3dr=zeros(npts(k),3); % the predicted 3D position with optimization.
    pltr3dr=zeros(npts(k),3); % the predicted 3D position with optimization.

    is_kv=zeros(npts(k),1);
    is_do=zeros(npts(k),1);
    pltis_do=zeros(npts(k),1);
    gAng=zeros(npts(k),1);
    chpt=zeros(npts(k),3);
    perror3dx=zeros(npts(k),1); % the prediction error in 3D
    pltperror3dx=zeros(npts(k),1); % the prediction error in 3D
    recon_error3d=zeros(npts(k),1); % the 3D error of estimation
    pltrecon_error3d=zeros(npts(k),1); % the 3D error of estimation
    perror2d=zeros(npts(k),1);  % the beam's eye view 2D error after prediction
    pltperror2d=zeros(npts(k),1);  % the beam's eye view 2D error after prediction

    nactual=nactual+npts(k); % the number of data points read.
    n_segment=0;             % the number of trajectories
    jj=1;                     % index of data points 
    j0=1;                    % index of the start point of a trajectory
    
 % the following stores correlation information in the original 3D data
 % for the first 4 seconds, used to screen trajectories
    kx=zeros(npts(k),3);
    
    extentx=zeros(npts(k),3);
    centerx=zeros(npts(k),3);
    centerpoint=[0 0 0];
    range_o3d=[0 0 0];
    range_r3d=[0 0 0];
    is_plane=0;
    linearity=ones(npts(k),1);
    n_domin=1;
    k1=0;
    k2=0;
    k3=0;  
    ko1=0;
    ko2=0;
    ko3=0;  
    for ii=1:npts(k)
        ko3d(ii,1)=o3d(k,ii,1);
        ko3d(ii,2)=o3d(k,ii,2);
        ko3d(ii,3)=o3d(k,ii,3);
        ko3d(ii,4)=o3d(k,ii,4);
    end  
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 % cycle the whole data in a subfolder(corresponding to a fraction)
while jj<npts(k) && npts(k)>points_period
% while jj<2500
        r3d0=zeros(runlength,3);  % the initial noise-added "original" 3D data
        is_jump=0;   % index indicating an abnormal trajectory due to 
                     % excessive velocity, which is to be discarded. 
%         checkpoint=[0 0 0]; % the center position in the initial 4 seconds.
        
        % Initializing: add detection error to the points of the first 4
        % seconds, calculate the extent, correlation and linearity of the 
        % motion for the first 4 seconds, and detect whether there is an 
        % unrealistic velocity(if yes, is_jump=1)
        if j0>npts(k)-points_period
            break
        end
        [r3d0,centerpoint,is_plane,linearity0,n_domin,k1,k2,k3,ko1,ko2,ko3,range_o3d0,range_r3d0]=...
            initializing_new(seg_err,j0,runlength,stdx,ko3d,sid);
        is_do(j0:j0+runlength-1)=ones(runlength,1);
        extentx(j0+runlength-1,:)=range_r3d0(:);
        centerx(j0+runlength-1,:)=[centerpoint(1) centerpoint(2) centerpoint(3)];
        kx(j0+runlength-1,:)=[k1 k2 k3];
        
        chpt(j0+runlength,:)=centerpoint(:);
        chpt(j0:j0+runlength-1,:)=zeros(runlength,3);  
        max_init=max(o3d(k,j0:j0+runlength-1,1));
        min_init=min(o3d(k,j0:j0+runlength-1,1));
        
        % Screen for unqualified trajectories of the first 4 seconds, and
        % skip them.
        % The conditions include: 1)poor linearity; 2)Extraordinarily large
        % or small amplitude of motion; 3)excessive velocity; 4) when the 
        % time of the remaining points is less than duration.

        while (range_o3d0(1)<0.1 || range_o3d0(1) >50 || ...
            range_o3d0(2)<0.1 || range_o3d0(2) >21 || ...
            range_o3d0(3)<0.2 || range_o3d0(3) >30 || ...
            abs(o3d(k,j0,1)-max_init)<0.01 || abs(o3d(k,j0,1)-min_init)<0.01 || ...%            linearity0< thres_linearity | 
            is_jump >0.1) && ...             
            j0<npts(k)-points_period-1
        
            j0=j0+1;
            [r3d0,centerpoint,is_plane,linearity0,~,k1,k2,k3,ko1,ko2,ko3,range_o3d0,range_r3d0]=...
                initializing_new(seg_err,j0,runlength,stdx,ko3d,sid);
            is_do(j0:j0+runlength-1)=ones(runlength,1);
            kx(j0+runlength-1,:)=[k1 k2 k3];
            max_init=max(ko3d(j0:j0+runlength-1,1));
            min_init=min(ko3d(j0:j0+runlength-1,1));
            extentx(j0+runlength-1,:)=range_r3d0(:);
            centerx(j0+runlength-1,:)=[centerpoint(1) centerpoint(2) centerpoint(3)];
            
            chpt(j0+runlength,:)=centerpoint(:);
            chpt(j0:j0+runlength-1,:)=zeros(runlength,3);           
        end 

        % if the time of the remaining points is less than duration, stop
        % the cycle.        

        
        % for the first 4 seconds, the noise-added data is directly
        % assigned to the reconstructed data, because they are measured
        % directly.
        r3d(j0:j0+runlength-1,1:3)=r3d0(1:runlength,1:3); 
        
        jj=runlength+j0;   % index of the first point to be reconstructed.
        is_discon=0;      % flag for interrupting the tracking.
        is_success=0;     % flag for complete tracking of a whole trajectory.
        is_drift=0;
        
        % initializing the projection points after the first 4 seconds.
        while is_discon<1                        
            % projecting and reconstructing the jth point.
            tic();
            [gAng(jj),is_do(jj)]=assign_angle_fixed(jj,j0,ko3d,angleset,onofftime);                       
            if is_do(jj)==1    
                [r3dx,is_discon,centerpoint,kx(jj-1,:),ko(jj-1,:),range_o3d,is_jump,linearity(jj-1)]= ...
                    proj_backproj_fixed(seg_err,ko3d,r3d,jj,j0,gAng(jj),runlength,sid,stdx);                
                extentx(jj-1,:)=range_o3d(:);
                centerx(jj-1,:)=[centerpoint(1) centerpoint(2) centerpoint(3)];
                % accepting the returned values.
                r3d(jj,1)=r3dx(1);
                r3d(jj,2)=r3dx(2);
                r3d(jj,3)=r3dx(3);           
            else
                r3d(jj,1)=NaN;
                r3d(jj,2)=NaN;
                r3d(jj,3)=NaN;
            end
            xx1=toc();
            timex(1)=timex(1)+xx1;
            if xx1>0.020
                fprintf('extra long estimation time:%f \n', xx1);
            end
            jj=jj+1;
            % if one segment is complete, move to next trajectory
            if jj>j0+points_period 
                is_discon=2;  % flag for terminating the cycle.
                is_success=1; % indicating a completed tracking
            end
        end % discon>1-... 

    % if the cycle of reconstruction is completed, make prediction and 
    % calculate the error & plot the graph
         if is_success>0.5 && jj<npts(k) 
             % the number of cycle steps up
             n_segment=n_segment+1;              
             % the calculation is done starting from the j0+n_lag, the
             % first n_lag points are neglected.
             for ik=1:n_lag
                 r3dr(j0+ik-1,1)=NaN;
                 r3dr(j0+ik-1,2)=NaN;
                 r3dr(j0+ik-1,3)=NaN;
                 r3dp(j0+ik-1,1)=NaN;
                 r3dp(j0+ik-1,2)=NaN;
                 r3dp(j0+ik-1,3)=NaN;
             end
% calculate the 3D extent of the trajectory in the first 4 seconds, for later reference.
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
             % start predicting
             is_used=[0 0 0];
             maxi=[0 0 0];
             %initializing the avarage
             av3d=NaN(points_period,3);
             for jx=j0+n_lag:jj-1        
%                  [gAng(jx-n_lag),is_do(jx-n_lag)]=assign_angle(k,jx-n_lag,j0,o3d,angleset,onofftime);    
                 if is_do(jx-n_lag)==1
                     av3d(jx-j0+1-n_lag,1)=r3d(jx-n_lag,1);
                     av3d(jx-j0+1-n_lag,2)=r3d(jx-n_lag,2);
                     av3d(jx-j0+1-n_lag,3)=r3d(jx-n_lag,3);                      
                     if jx>j0+n_lag+1
                         if is_do(jx-n_lag-2)==1
                             av3d(jx-j0-n_lag,1)=(r3d(jx-n_lag-2,1)+2*r3d(jx-n_lag-1,1)+r3d(jx-n_lag,1))/4;
                             av3d(jx-j0-n_lag,2)=(r3d(jx-n_lag-2,2)+2*r3d(jx-n_lag-1,2)+r3d(jx-n_lag,2))/4;
                             av3d(jx-j0-n_lag,3)=(r3d(jx-n_lag-2,3)+2*r3d(jx-n_lag-1,3)+r3d(jx-n_lag,3))/4;                
                         elseif is_do(jx-n_lag-1)==1
                             av3d(jx-j0-n_lag,1)=r3d(jx-n_lag-1,1);
                             av3d(jx-j0-n_lag,2)=r3d(jx-n_lag-1,2);
                             av3d(jx-j0-n_lag,3)=r3d(jx-n_lag-1,3);
                         end
                     end 
                 end

                 if jx>j0+runlength && n_lag>0.1 && is_do(jx)==1
                     % checkpoint updated
%                      chpt(jx,:)=chpt(jx-1,:);
                     % calculate the immediate past range
                     extentx=near_range3d(jx,j0,runlength,n_lag,r3d);                     
                     % calculate the difference between the point n_log
                     % before and the current point, for prediction. Note the
                     % input information is from the points at least n_lag
                     % before. 
                     
                     %% for the first two points in the treatment interval, use
                     %% the last two points estimated as the predicted
                     %% position
                     if is_do(jx-2)==0 || is_do(jx-1)==0 || is_do(jx-n_lag)==0                        
                         r3dr(jx,:)=[r3d(jx-points_transition-3,1) r3d(jx-points_transition-3,2) r3d(jx-points_transition-3,3)];
                     elseif is_do(jx-n_lag-1)~=0
                         tic();
                         [incr_3d(jx-n_lag,1) incr_3d(jx-n_lag,2) incr_3d(jx-n_lag,3)]=predicting_fixed(is_pre,factor1,factor2,n_lag,jx-n_lag,j0,extentx,av3d);                     
                         % Add the difference to the position n_lag points before 
                         % to get the predicted position of the current point.
                         % If the calculated difference is too large, replace
                         % it with a threshold value, which is based on the
                         % amplitude of motion in the first 4 seconds. 

                         [r3dp(jx,1) r3dp(jx,2) r3dp(jx,3)]=prediction_screening_fixed(jx,n_lag,r3d,incr_3d,extentx,factor3);
                         % optimize the prediction result with the updated 
                         % correlation information n_lag points before 
                         if jx>j0+n_lag+runlength
                             r3dr(jx,:)=prediction_optimizing_new(jx,r3dp,stdx,centerx(jx-n_lag,:),kx(jx-n_lag,:),linearity(jx-n_lag)); % true realtime  
                         else
                             r3dr(jx,:)=[r3dp(jx,1) r3dp(jx,2) r3dp(jx,3)];
                         end
                         xyy=toc();
                         if xyy>0.030
                             fprintf('extraordinary prediction time =%f ! \n', xyy);
                         end                         
                     else
                         r3dr(jx,1:3)=r3d(jx-n_lag,1:3);
                     end 
                     [recon_error3d(jx) perror3dx(jx) perror2d(jx)]=error_calculation(k,jx,gAng(jx),sid,stdx,o3d,r3d,r3dr);
                      npoints(k)=npoints(k)+1;   
%                       if isnan(recon_error3d(jx)) | isnan(perror3dx(jx)) | isnan(perror2d(jx))
%                           pause(0.1);
%                       end
                 elseif is_do(jx)==1 % of jx>j0+runlength ...
                      % if the point position is measured in the first 4 seconds,
                      % no prediction is made, use the latest point estimated  
                     r3dr(jx,1:3)=r3d(jx-n_lag,1:3);
                     % calculate the projection of estimated position.
                     [recon_error3d(jx) perror3dx(jx) perror2d(jx)]=error_calculation(k,jx,gAng(jx),sid,stdx,o3d,r3d,r3dr);
                     npoints(k)=npoints(k)+1;
                 else
                     r3dr(jx,1:3)=[NaN NaN NaN];
                     recon_error3d(jx)=NaN;
                     perror3dx(jx)=NaN;
                     perror2d(jx)=NaN;
                 end % of jx>j0+runlength ...

                 % the number of points predicted increases by one.
                 % 
                 if is_do(jx)~=0
                     aberror2(k)=aberror2(k)+perror3dx(jx)*perror3dx(jx);
                     ab_dev2(k)=ab_dev2(k)+perror2d(jx)*perror2d(jx);
                 end
             end % of "for jx=j0+n_lag:jj-1"
               
             % calculate the variation of the trajectory itself.
             od=zeros(jj-j0,3);
             % construct a new array for calculating the standard deviation
             % of the trajectory.
             od(:,1)=ko3d(j0:jj-1,1);
             od(:,2)=ko3d(j0:jj-1,2);
             od(:,3)=ko3d(j0:jj-1,3);
             % sum up the square of all standard deviation.
             std_deviation=sqrt(std(od(1:jj-j0,1))*std(od(1:jj-j0,1))+...
                 std(od(1:jj-j0,2))*std(od(1:jj-j0,2))+...
                 std(od(1:jj-j0,3))*std(od(1:jj-j0,3)));
             std_dev2(k)=std_deviation*std_deviation*(jj-j0-1)+std_dev2(k); 
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % draw the figure of the error and predicated vs. actual
             % position.             
             xcoord=zeros(points_period,1);
             
             std1=sqrt(ab_dev2(k)/(points_period-points_transition*(length(angleset)-1))); % average 2D error with prediction
                         
             percent95_3dx=get_percentile_error_fixed(perror3dx,0.95,j0+n_lag,j0+points_period-1);
             percent95_2dx=get_percentile_error_fixed(perror2d,0.95,j0+n_lag,j0+points_period-1);
             percent99_3dx=get_percentile_error_fixed(perror3dx,0.99,j0+n_lag,j0+points_period-1);
             percent99_2dx=get_percentile_error_fixed(perror2d,0.99,j0+n_lag,j0+points_period-1);
             percent95_3d(k)=percent95_3d(k)+percent95_3dx;
             percent95_2d(k)=percent95_2d(k)+percent95_2dx;
             percent99_3d(k)=percent99_3d(k)+percent99_3dx;
             percent99_2d(k)=percent99_2d(k)+percent99_2dx;
             if percent95_2dx>3.0
                 count95_3=count95_3+1;
                 if percent95_2dx>4.0
                     count95_4=count95_4+1;
                 end                 
             end
             onoff=ones(points_period,1);
             kxi=1;
             for ii=1:points_period
                 xcoord(ii)=ko3d(ii+j0-1,4)-ko3d(j0,4);                     
                 if abs(xcoord(ii)-onofftime(kxi))<0.0385*incrt/2                        
                     onoff(kxi)=ii;  
                     if kxi~=length(onofftime)
                         kxi=kxi+1;
                     end                     
                 end                 
             end
             indxx=[onoff(1):onoff(2), onoff(2):8:onoff(3), ...
                 onoff(3):onoff(4), onoff(4):8:onoff(5), ...
                 onoff(5):onoff(6), onoff(6):8:onoff(7), ...
                 onoff(7):onoff(8), onoff(8):8:onoff(9), ...
                 onoff(9):onoff(10)];  
             kj=1;
             while kj<length(indxx)
                 plto3d(kj,1)=o3d(k,j0-1+indxx(kj),1);
                 plto3d(kj,2)=o3d(k,j0-1+indxx(kj),2);
                 plto3d(kj,3)=o3d(k,j0-1+indxx(kj),3);
                 pltr3d(kj,1)=r3d(j0-1+indxx(kj),1);
                 pltr3d(kj,2)=r3d(j0-1+indxx(kj),2);
                 pltr3d(kj,3)=r3d(j0-1+indxx(kj),3);
                 pltr3dr(kj,1)=r3dr(j0-1+indxx(kj),1);
                 pltr3dr(kj,2)=r3dr(j0-1+indxx(kj),2);
                 pltr3dr(kj,3)=r3dr(j0-1+indxx(kj),3);
                 pltperror3dx(kj)=perror3dx(j0-1+indxx(kj));
                 pltperror2d(kj)=perror2d(j0-1+indxx(kj));
                 pltrecon_error3d(kj)=recon_error3d(j0-1+indxx(kj));
                 pltis_do(kj)=is_do(j0-1+indxx(kj));
                 kj=kj+1;
             end

            conditionx=std1<1.47 & std_deviation>2.9; %& k==1 & j0 > 11750 & j0<11800;
            modex=2;            
%             if num_plot<0.5
%                 num_plot=plotfigures_fixed(conditionx,modex,n_segment,k,j0,indxx,...
%                     plto3d,pltr3dr,pltr3d,pltperror2d,pltperror3dx,pltrecon_error3d,std1,dbfolders);
%                 n_segment=n_segment+200;
%                 modex=1;
%                 num_plot=plotfigures_fixed(conditionx,modex,n_segment,k,j0,indxx,...
%                     plto3d,pltr3dr,pltr3d,pltperror2d,pltperror3dx,pltrecon_error3d,std1,dbfolders);
%                 n_segment=n_segment-200;
%             end                 
         end % of "is_success>0.5 & jj<npts(k)."
         is_success=0;  % reset the flag of completing a trajectory      
         j0=jj; % one cycle is finished, shift the starting point      
    end % of while jj<npts(k)
    if npts(k)>points_period
        n_seg=n_seg + n_segment; % the number of trajectories increases by n_segment.
        if n_segment~=0
            percent95_3d(k)=percent95_3d(k)/n_segment;
            percent95_2d(k)=percent95_2d(k)/n_segment;
            percent99_3d(k)=percent99_3d(k)/n_segment;
            percent99_2d(k)=percent99_2d(k)/n_segment;
            std_dev2(k)=std_dev2(k)/npoints(k);
            aberror2(k)=aberror2(k)/npoints(k);
            ab_dev2(k)=ab_dev2(k)/npoints(k);     
        else
            percent95_3d(k)=0;
            percent95_2d(k)=0;
            percent99_3d(k)=0;
            percent99_2d(k)=0;
            std_dev2(k)=0;
            aberror2(k)=0;
            ab_dev2(k)=0;           
        end
    end
end % of k...
std_dev=sum(sqrt(std_dev2(1:ndir)))/ndir; % total standard deviation of the patient.
rms=sum(sqrt(aberror2(1:ndir)))/ndir;     % total rms 3D error
devp=sum(sqrt(ab_dev2(1:ndir)))/ndir;     % total 2D error
pct95_3d=mean(percent95_3d(1:ndir));
pct95_2d=mean(percent95_2d(1:ndir));
pct99_3d=mean(percent99_3d(1:ndir));
pct99_2d=mean(percent99_2d(1:ndir));

% indicate the analysis of one patient is finished
fprintf('the number of trajectories = %f \n', n_seg); 
fprintf('the number of fractions = %f \n', ndir); 
return