function parameter_control()
% this file is to the final execution file, which could compare different
% parameters.
clear all
close all
% is_pre = 1 if prediction is applied, 0 otherwise.
is_pre = 1; % conditions(1);
% is_block = 1 if block is simulated.
is_block = 0; % conditions(2);
% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. selected range = full extent/factor1
% in that direction in the first 4 seconds.
% % factor1 = 4.0; %conditions(3);
% % 
% % % factor2 specifies the range of velocity within which the points are
% % % selected for multivariate regression. selected range = max
% % % value*factor2
% % % in that direction in the first 4 seconds.
% % factor2 = 0.20; %conditions(4);
factor1 = 4.0;
factor2 = 0.2;

r2threshold = 0.10;
num_train = 12;
% factor3 specifies the ratio of the threshold controlling the predicted 
% changes of positions to the extent in their respective direction.
factor3 = 1.4; % conditions(5);
seg_err = 0;  %conditions(6); % specifies the rms error of detection(=0.5mm).
lag_time = 310; %conditions(7); % the latent time in milliseconds
% lag_time = 0; %conditions(7); % the latent time in milliseconds
thres_linearity = 0.01;% conditions(8); % the threshold of correlation coefficient for the motion of 
% the first 4 seconds. If the threshould is not reached, the data of the 4
% second is discarded. The program moves to the data of next 4 seconds.
amp = 4.0; % conditions(9); % specifies the maximum ratio of motion amplitude of 
                                % later data points to that of the first 4 seconds.
duration_ = 72; % conditions(10); % the time interval for one cycle of arc therapy (=72)
init_gAng = -179; % conditions(11); % initial angle(=-179)
end_gAng = 179; % conditions(12); % final angle(=+179)
% runlength: the number of initial points whose actual 3D coordinates
%  were measured in the first few seconds(~ 26)
init_time = 4; % conditions(13);
incrt = 2; % conditions(14);
runlength = floor(init_time/(incrt*0.0385))+1;
% one sample point are taken ever incrt data points, whose 
% period is approximately 0.038545s.incrt(= 4)
sid = 1500; % conditions(15);  the Source-Imager distance.(=1500)
stdx = 1000;%  conditions(16);  the source-isocenter distance.(=1000)
min_SI = 0.2;% conditions(17); % the min SI amplitude for the first few seconds ~0.2;
max_SI = 60; % conditions(18); % the max SI amplitude for the first few seconds ~50;
min_AP = 0.1; % conditions(19); % the min AP amplitude for the first few seconds ~0.1;
max_AP = 30; % conditions(20); % the max AP amplitude for the first few seconds ~27;
min_LR = 0.1; % conditions(21); % the min LR amplitude for the first few seconds ~0.1;
max_LR = 30; % conditions(22); % the max LR amplitude for the first few seconds ~20;
marginx = 0.3;

conditions =[is_pre,is_block,factor1, factor2, factor3,seg_err,lag_time,...
    thres_linearity, amp, duration_,init_gAng,end_gAng, 4, incrt, sid, stdx,...
    min_SI,max_SI,min_AP,max_AP,min_LR,max_LR,r2threshold,num_train,marginx];
delete('e:\research\Liuwu\stat_datax2');
patient_selection(conditions);
return