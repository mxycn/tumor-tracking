% function [k_segment, std_k2,pre3d_k2,pro2d_k2,pct95_3d,pct95_2d,pct99_3d,pct99_2d] = ...
%     simulation_fraction(npts,ko3d,conditions)
function [k_segment, std_k2,pre3d_k2,pro2d_k2,threedf,twodf,parametersf,periodic_contrastf] = ...
    simulation_fraction(npts,ko3d,conditions)

% this is used to simulate the treatment with data of a fraction.

% n_segment: the number of different segments(trajectories) in one
% fraction.
% npoints: the number of points with prediction.
% aberror2: the average of the square of 3D error
% ab_dev2: the average of the square of 2D error
% std_dev2: the sum of the square of standard deviation, indicating the
% amplitude of motion

%% initializing the parameters and arrays
incrt = conditions(14);
duration_ = conditions(10); % the time interval for one cycle of arc therapy (=72)
init_gAng = conditions(11); % initial angle(=-179)
end_gAng = conditions(12); % final angle(=+179)
angular_speed=abs(end_gAng-init_gAng)/duration_; % the rate of gantry rotation 
ang_incrt=angular_speed*incrt*0.038545; % the angle interval between data points.
points_period=floor(358.0/ang_incrt);   % the number of points in a trajectory.        
init_time = conditions(13);
runlength = floor(init_time/(incrt*0.0385))+1;
lag_time = conditions(7);
n_lag=round(lag_time/(incrt*38.5));

% true range and detected range in the first few seconds and the whole range.
range_o3d0 = zeros(1,3);
% range_o3d = zeros(1,3);

thres_linearity = conditions(8); % the threshold of correlation coefficient for the motion of 
% the first 4 seconds. If the threshould is exceeded, the data of the 4
% second is discarded. The program moves to the data of next 4 seconds.

amp = conditions(9); % specifies the maximum ratio of motion amplitude of 
                                % later data points to that of the first 4 seconds.
min_SI = conditions(17); % the min LR amplitude for the first few seconds ~0.2;
max_SI = conditions(18); % the min LR amplitude for the first few seconds ~50;
min_AP = conditions(19); % the min AP amplitude for the first few seconds ~0.1;
max_AP = conditions(20); % the min LR amplitude for the first few seconds ~27;
min_LR = conditions(21); % the min LR amplitude for the first few seconds ~0.1;
max_LR = conditions(22); % the min LR amplitude for the first few seconds ~20;
                               
std_k2 = 0;
pre3d_k2 = 0;
pro2d_k2 = 0;
k_segment = 0;
pct95_3d = 0;
pct95_2d = 0;
pct99_3d = 0;
pct99_2d = 0;
j0=1; % the index of starting points of each segments 
%% simulating ...
while j0<= npts-points_period 
    % test whether the coming points_period points meet the condition
    is_success = 0;
%     tic();
    while is_success == 0 && j0<= npts-points_period 
        xx0=range(ko3d(j0:j0+runlength-1,1:3),1);
        xx=range(ko3d(j0:j0+points_period-1,1:3),1);
        range_o3d0=xx0(1,:);
        
        % range_o3d=xx(1,:);
        % the ratio of range between the whole trajectory and that of the
        % first few seconds.
        ampx = max(xx./xx0);
        
        % To calculate the range of the speed
        maxv = max(max(abs(ko3d(j0+1:j0+points_period-1,1:3)-ko3d(j0:j0+points_period-2,1:3)),[],1));
        time_interval = max(abs(ko3d(j0+1:j0+points_period-1,4)-ko3d(j0:j0+points_period-2,4)));        
        % the calculate the linearity of the whole trajectory.
        [coeffo,scoreo,rootso] = princomp(ko3d(j0:j0+points_period-1,1:3));        
        linearity = rootso(1)/sum(rootso);
        
        % find the maximum and minimum SI cooridinate for the first few
        % seconds.
        max_init=max(ko3d(j0:j0+runlength-1,1));
        min_init=min(ko3d(j0:j0+runlength-1,1)); 
        % The conditions include: 1)poor linearity; 2)Extraordinarily large
        % or small amplitude of motion; 3)excessive velocity; 4) when the 
        % time of the remaining points is less than duration; 5) initial point
        % is on the ends of the trajectory.
        if time_interval > 0.19 ||  maxv > 6.5
            j0 = j0 + runlength;
        elseif (range_o3d0(1)<min_SI || range_o3d0(1) >max_SI || ...
            range_o3d0(2)<min_AP || range_o3d0(2) >max_AP || ...
            range_o3d0(3)<min_LR || range_o3d0(3) >max_LR ) 
            j0 = j0 + 3;
        elseif (abs(ko3d(j0,1)-max_init) < 0.01 || abs(ko3d(j0,1)-min_init)< 0.01 || ...
                ampx > amp || linearity < thres_linearity)
            j0 = j0 + 1;
        else
            is_success = 1;
            j0 = j0 + 1;
        end
    end % of "is_success == 0"
    if j0 > npts-points_period + 1
       break
    end
    j0 = j0-1;
    % data localization 
    so3d = zeros(points_period,3);
    so3d(:,:)=ko3d(j0:j0+points_period-1,1:3);
%     xt1=toc();
    % simulate the a segment, i.e., a trajectory
%     tic();
    [std_dev2,aberror2,ab_dev2,parameters,threed,twod,WaperR_adjusted] = simulation_segment(so3d,conditions);    
    prct95_3d = 0;
    prct95_2d = 0;
    prct99_3d = 0;
    prct99_2d = 0;
%     [std_dev2,aberror2,ab_dev2,prct95_3d,prct95_2d,prct99_3d,prct99_2d] = simulation_segment(so3d,conditions);    
%     xt2=toc();
    k_segment = k_segment + 1
    if k_segment==1
        threedf = threed;
        twodf = twod;
        parametersf = parameters;
        periodic_contrastf=[std_dev2,aberror2,ab_dev2,WaperR_adjusted];
    else
        threedf = [threedf;threed];
        twodf = [twodf;twod];
        parametersf = [parametersf; parameters];
        periodic_contrastf=[periodic_contrastf;std_dev2,aberror2,ab_dev2,WaperR_adjusted];
    end
    j0 = j0 + points_period;
%     fprintf('the number of segments = %u, time1 = %f vs. time2 = %f \n', k_segment, xt1, xt2); 
%     fprintf('standard deviation = %f, 3D error = %f vs. 2D error %f \n', sqrt(std_dev2/points_period),sqrt(aberror2/(points_period-n_lag)),sqrt(ab_dev2/(points_period-n_lag))); 
    std_k2 = std_k2 + std_dev2;
    pre3d_k2 = pre3d_k2 + aberror2;
    pro2d_k2 = pro2d_k2 + ab_dev2;
    
    pct95_3d = pct95_3d + prct95_3d;
    pct95_2d = pct95_2d + prct95_2d;
    pct99_3d = pct99_3d + prct99_3d;
    pct99_2d = pct99_2d + prct99_2d;
end 
%% convert the errors to their standard square form
std_k2 = std_k2/k_segment;
pre3d_k2 = pre3d_k2/k_segment;
pro2d_k2 = pro2d_k2/k_segment;    
pct95_3d = pct95_3d/k_segment; 
pct95_2d = pct95_2d/k_segment; 
pct99_3d = pct99_3d/k_segment; 
pct99_2d = pct99_2d/k_segment; 
% pause(1);

