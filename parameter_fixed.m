function parameter_fixed()
% this file is to the final execution file, which could compare different
% parameters.

% is_pre = 1 if prediction is applied, 0 otherwise.

% factor1 specifies the range of coordinate within which the points are
% selected for multivariate regression. factor1=selected range/full extent
% in that direction in the first 4 seconds.

% factor2 specifies the range of velocity within which the points are
% selected for multivariate regression. factor2=selected range/full extent
% in that direction in the first 4 seconds.

% factor3 specifies the ratio of the threshold controlling the predicted 
% changes of positions to the extent in their respective direction.

% n_lag is the number of latent points, n=3 corresponds to 460ms, n=2
% cooresponds to 310ms.

% init_length(=26, corresponding to 4 seconds) is the number of points for 
% initial full 3D reconstruction. 

% linearity presents the threshold of correlation coefficient for the motion of 
% the first 4 seconds. If the threshould is surpassed, the data of the 4
% second is discarded. The program moves to the data of next 4 seconds.

% amp specifies the maximum ratio of motion amplitude of later data points
% to that of the first 4 seconds.
% patient_selection(incrt,is_pre,2.7,16,2.5,seg_err,lag_time,init_time,linearity,multiple); % with prediction and with seg.
angleset=[35 110 180 250 325]-180;
onofftime=[0 20 70 90 140 160 210 230 280 300];
num_plot=0;
% while num_plot<0.5
    num_plot=patient_selection_fixed(4,1,0,2.6,16,2.6,0.5,460,2.0,0.01,5,angleset,onofftime); % with prediction
% end
% patient_selection_fixed(4,1,0,2.6,16,3.3,0.5,310,4.0,0.01,5,angleset,onofftime); % with prediction
% patient_selection_fixed(4,0,0,2.6,16,2.6,0.5,460,4.0,0.01,5,angleset,onofftime); % with prediction and with seg.
% patient_selection_fixed(4,0,0,2.6,16,3.3,0.5,310,4.0,0.01,5,angleset,onofftime); % with prediction

% patient_selection_fixed(4,0,0,2.6,16,2.6,0.5,0.0,4.0,0.01,5,angleset,onofftime); % with prediction and with seg.
% patient_selection_fixed(4,0,0,2.6,16,3.3,0.5,310,4.0,0.01,5,angleset,onofftime); % with prediction
% patient_selection_fixed(4,0,0,2.6,16,3.3,0.5,0,4,0.01,5,angleset,onofftime); % with prediction and with seg.

% patient_selection_fixed(4,1,0,2.6,16,2.7,0.5,460,4,0.01,5); % with prediction and with

% seg.

return