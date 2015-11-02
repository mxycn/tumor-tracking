function [original,npts]=analyze_datax()
%
% analyze_data 
% read stat_datax, which contains all errors data and their matching
% parameters. And sort the errors according to their parameters which
% themselves have been binned into descrete numbers.
%
% OUTPUT: the averages and standard deviations for different bins of parameters.
% In the first step the parameters are the square of cosine and the cost
% function. The are both divided into 10 bins.
% Programmed by: Huagang Yan, yanhg@ccmu.edu.cn
% Aug 10, 2010
%--------------------------------------------------------------------------

% the upperlimit of npoints is 400,000,

original=zeros(400000,6);% since each folder may contain different numbers of data points
                           % here we set an upper limit.
sorted=zeros(10,10,40000);
aves=zeros(10,10);
stds=zeros(10,10);
avex=zeros(100,1);
stdx=zeros(100,1);

inds=ones(10,10);                         
%% read data file
fid = fopen('e:\research\Liuwu\stat_datax2','r');
ki=1;
while ~feof(fid) 
    tline = fgetl(fid); % after reading each line, move onto the next line.
        if ~ischar(tline), break, end
        [xT yT zT aT bT cT] = strread(tline,'%f %f %f %f %f %f',1);
        original(ki,:) = [sqrt(xT) sqrt(yT) zT aT bT cT];
        ki=ki+1;                  
end % of while statement
fclose(fid);
% read data file done...
%% group the data into 100 bins
% create discrete parameter values
cos2 = zeros(10,1);
cost = zeros(11,1);
cos2x = zeros(100,1);
costx = zeros(100,1);
cos2x2 = zeros(100,1);
preerr = zeros(11,1);
preerrx = zeros(100,1);
% divide each parameter into 10 bins.
for i = 1:10
    cos2(i) = (i-1)*0.1;
    preerr(i) = (i-1)*0.2;
end
preerr(11) = 2.5;
% classify the data according to the square of cosine theta(where theta is the angle 
% between the principal axis and treatment beam, and the cost function.
% The cost function may exceed 1, in that case, the data points will not be
% considered.
for k = 1:ki-1
    for i = 1:9
        for j = 1:10
            if original(k,2)>cos2(i) && original(k,2)<= cos2(i+1) && original(k,3)>preerr(j) && original(k,3)<=preerr(j+1) 
                sorted(i,j,inds(i,j)) = original(k,6);
                inds(i,j)=inds(i,j)+1;
            end
        end
    end
    % the cosine function has a limit. 
    for j = 1:10
        if original(k,2)>cos2(10) && original(k,3)>preerr(j) && original(k,3)<=preerr(j+1) 
            sorted(10,j,inds(10,j)) = original(k,6);
            inds(10,j)=inds(10,j)+1;
        end
    end
end
for i = 1:10
    for j = 1:10        
        stds(i,j)=std(sorted(i,j,1:inds(i,j)));
        aves(i,j)=mean(sorted(i,j,1:inds(i,j)));
    end
end
% flatten the stds and aves into one dimensional array, as well as the parameters into 1D array.
for i = 1:100
    cos2x(i) = cos2(ceil(i/10)); % the square of cosine theta, where theta is the angle between the principal axis and treatment bean.   
    cos2x2(i) = cos2x(i)^2;
    if mod(i,10)>0 
        j = mod(i,10);
        preerrx(i) = preerr(j);        
    else
        j = 10;
        preerrx(i) = preerr(10);        
    end
    avex(i) = aves(ceil(i/10.0),j);
    stdx(i) = stds(ceil(i/10.0),j);        
end
% [regr1 r1]=two_var_regress(100,cos2x(1:100),costx(1:100),avex(1:100));
[bb1,bint1,r1,rint1,stats1] = regress(avex,[ones(100,1) cos2x preerrx cos2x2] );
% [regr2 r2]=two_var_regress(100,cos2x(1:100),costx(1:100),stdx(1:100));
[bb2,bint2,r2,rint2,stats2] = regress(stdx,[ones(100,1) cos2x preerrx cos2x2] );

[bb3,bint3,r3,rint3,stats3] = regress(avex,[ones(100,1) cos2x preerrx] );
% [regr2 r2]=two_var_regress(100,cos2x(1:100),costx(1:100),stdx(1:100));
[bb4,bint4,r4,rint4,stats4] = regress(stdx,[ones(100,1) cos2x preerrx] );
pause(1);
return
          


