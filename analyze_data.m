function [original,npts]=analyze_data()
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

original=zeros(400000,9);% since each folder may contain different numbers of data points
                           % here we set an upper limit.
sorted1=zeros(10,10,30000);
sorted2=zeros(10,10,30000);
sorted3=zeros(10,10,30000);
sorted4=zeros(10,10,30000);
sorted5=zeros(10,10,30000);
sorted6=zeros(10,10,30000);
stds1=zeros(10,10);
stds2=zeros(10,10);
stds3=zeros(10,10);
stds4=zeros(10,10);
stds5=zeros(10,10);
stds6=zeros(10,10);
avex=zeros(100,1);
stdx1=zeros(100,1);
stdx2=zeros(100,1);
stdx3=zeros(100,1);
stdx4=zeros(100,1);
stdx5=zeros(100,1);
stdx6=zeros(100,1);

inds1=ones(10,10);                   
inds2=ones(10,10); 
%% read data file
fid = fopen('e:\research\Liuwu\stat_datax3','r');
ki=1;
while ~feof(fid) 
    tline = fgetl(fid); % after reading each line, move onto the next line.
        if ~ischar(tline), break, end
        [v cos2 pree preq costx aT bT cT xT yT zT] = strread(tline,'%f %f %f %f %f %f %f %f %f %f %f',1);
        original(ki,:) = [sqrt(cos2) sqrt(costx) sqrt(pree) aT bT cT xT yT zT];
        ki=ki+1;                  
end % of while statement
fclose(fid);
% read data file done...
%% group the data into 100 bins
% create discrete parameter values
cos2 = zeros(10,1);
cost = zeros(11,1);
pre = zeros(11,1);
cos2x = zeros(100,1);
cos2x2 = zeros(100,1);
costx = zeros(100,1);
costx2 = zeros(100,1);
prex = zeros(100,1);
prex2 = zeros(100,1);
% divide each parameter into 10 bins.
for i = 1:10
    cos2(i) = (i-1)*0.1;
end
for i = 1:11
    cost(i) = (i-1)*0.1;
end
for i = 1:11
    pre(i) = (i-1)*0.20;
end
% classify the data according to the square of cosine theta(where theta is the angle 
% between the principal axis and treatment beam, and the cost function.
% The cost function may exceed 1, in that case, the data points will not be
% considered.
for k = 1:ki-1
    for i = 1:9
        for j = 1:10
            if original(k,1)>cos2(i) && original(k,1)<= cos2(i+1) && original(k,2)>cost(j) && original(k,2)<=cost(j+1) 
                sorted1(i,j,inds1(i,j)) = original(k,4);
                sorted2(i,j,inds1(i,j)) = original(k,5);
                sorted3(i,j,inds1(i,j)) = original(k,6);
                inds1(i,j)=inds1(i,j)+1;
            end
            if original(k,1)>cos2(i) && original(k,1)<= cos2(i+1) && original(k,3)>pre(j) && original(k,3)<=pre(j+1) 
                sorted4(i,j,inds2(i,j)) = original(k,7);
                sorted5(i,j,inds2(i,j)) = original(k,8);
                sorted6(i,j,inds2(i,j)) = original(k,9);
                inds2(i,j)=inds2(i,j)+1;
            end
        end
    end
    % the cosine function has a limit. 
    for j = 1:10
        if original(k,1)>cos2(10) && original(k,2)>cost(j) && original(k,2)<=cost(j+1) 
            sorted1(10,j,inds1(10,j)) = original(k,4);
            sorted2(10,j,inds1(10,j)) = original(k,5);
            sorted3(10,j,inds1(10,j)) = original(k,6);
            inds1(10,j)=inds1(10,j)+1;
        end
        if original(k,1)>cos2(10) && original(k,3)>pre(j) && original(k,3)<=pre(j+1) 
            sorted4(10,j,inds2(10,j)) = original(k,7);
            sorted5(10,j,inds2(10,j)) = original(k,8);
            sorted6(10,j,inds2(10,j)) = original(k,9);
            inds2(10,j)=inds2(10,j)+1;
        end
    end
end
for i = 1:10
    for j = 1:10        
        stds1(i,j)=std(sorted1(i,j,1:inds1(i,j)));
        stds2(i,j)=std(sorted2(i,j,1:inds1(i,j)));
        stds3(i,j)=std(sorted3(i,j,1:inds1(i,j)));
        stds4(i,j)=std(sorted4(i,j,1:inds2(i,j)));
        stds5(i,j)=std(sorted5(i,j,1:inds2(i,j)));
        stds6(i,j)=std(sorted6(i,j,1:inds2(i,j)));
%         aves(i,j)=mean(sorted(i,j,1:inds(i,j)));
    end
end
% flatten the stds and aves into one dimensional array, as well as the parameters into 1D array.
for i = 1:100
    cos2x(i) = cos2(ceil(i/10))+0.05; % the square of cosine theta, where theta is the angle between the principal axis and treatment bean.   
    cos2x2(i) = cos2x(i)^2;
    if mod(i,10)>0 
        j = mod(i,10);
        costx(i) = cost(j)+0.05;        
        prex(i) = pre(j)+0.1;
        costx2(i)=costx(i)^2;
        prex2(i)=prex(i)^2;
    else
        j = 10;
        costx(i) = cost(10)+0.05;  
        prex(i) = pre(10)+0.1;
        costx2(i)=costx(i)^2;
        prex2(i)=prex(i)^2;
    end
%     avex(i) = aves(ceil(i/10.0),j);
    stdx1(i) = stds1(ceil(i/10.0),j);        
    stdx2(i) = stds2(ceil(i/10.0),j); 
    stdx3(i) = stds3(ceil(i/10.0),j); 
    stdx4(i) = stds4(ceil(i/10.0),j);        
    stdx5(i) = stds5(ceil(i/10.0),j); 
    stdx6(i) = stds6(ceil(i/10.0),j);
end
% for i = 1:150
%     cos2x(i) = cos2(ceil(i/15)); % the square of cosine theta, where theta is the angle between the principal axis and treatment bean.   
%     cos2x2(i) = cos2x(i)^2;
%     if mod(i,15)>0 
%         j = mod(i,15);
%         costx(i) = cost(j);        
%     else
%         j = 15;
%         costx(i) = cost(15);        
%     end
%     stdx4(i) = stds4(ceil(i/10.0),j);        
%     stdx5(i) = stds5(ceil(i/10.0),j); 
%     stdx6(i) = stds6(ceil(i/10.0),j);       
% end

subplot(2,3,1)
mesh(cos2(1:10),cost(1:10),stds1);
subplot(2,3,2)
mesh(cos2(1:10),cost(1:10),stds2);
subplot(2,3,3)
mesh(cos2(1:10),cost(1:10),stds3);
subplot(2,3,4)
mesh(cos2(1:10),cost(1:10),stds4);
subplot(2,3,5)
mesh(cos2(1:10),cost(1:10),stds5);
subplot(2,3,6)
mesh(cos2(1:10),cost(1:10),stds6);
% [bb2,bint2,r2,rint2,stats1] = regress(stdx2,[ones(100,1) cos2x cos2x2 costx] );
% [bb3,bint3,r3,rint3,stats3] = regress(stdx3,[ones(100,1) cos2x cos2x2 costx] );
% [bb4,bint4,r4,rint4,stats4] = regress(stdx4,[ones(100,1) cos2x cos2x2 prex] );
% [bb5,bint5,r5,rint5,stats5] = regress(stdx5,[ones(100,1) cos2x cos2x2 prex] );
[bb5,bint5,r5,rint5,stats5] = regress(stdx6,[ones(100,1) cos2x prex costx] );
[bb6,bint6,r6,rint6,stats6] = regress(stdx6,[ones(100,1) cos2x cos2x2 prex] );
pause(1);
return
          


