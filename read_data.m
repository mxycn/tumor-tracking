function [o_3D,npts]=read_data(ndir,incrt,path_in)
%
% read_data 
% read Modeler.log, which contains target positions estimated
%    by Synchrony Model {time,xT,yT,zT} at 25 Hz.
%    Note that it was processed from the original Modeler.log file.
%    The details is described in the PMB paper [2008;53:3623-3640].
%
% INPUT: the directory containing the data, e.g. 'D:\research\LiuWu\database\DB01'
% OUTPUT: the number of the points, an array of t_data(starting from zero),
% an 3D arrry of coordinates of the target.
% Programmed by: Huagang Yan, yanhg@ccmu.edu.cn
% Aug 10, 2010
%--------------------------------------------------------------------------

% the upperlimit of npoints is 300,000,

o_3D=zeros(ndir,300000,4); % since each folder may contain different numbers of data points
                           % here we set an upper limit.
% Define constants
% Need to set the directory

ModelerDataDir = path_in;
cd(ModelerDataDir); 
dbfolders=dir('DB*');

% the number of data points in each subfolder
npts=zeros(ndir,1);
% ********************************************
% ***** Loop over all patient data files *****
% ********************************************

% for folderIndex = 1:1
for folderIndex = 1:ndir
    cd(ModelerDataDir); 
    cd(dbfolders(folderIndex).name);
    subfolder = dir('DB*');
    cd(subfolder.name);
    PathName = cd;
    
    % ***************************
    % 1) read ModelPoints.log file
    % ***************************
    fprintf('%s\n',dbfolders(folderIndex).name);
    fid = fopen('ModelerExt.log','r');
    dataIndex = 1;
    while ~feof(fid)
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        dataIndex = dataIndex + 1;
    end % of while

    nPoints = dataIndex-2; %减2为了排除降采样后采到最末一点     
    npts(folderIndex)=floor(nPoints/incrt);%%%%%%%%%%%%%%%% second output
    
%    nData = zeros(nModelePoints,1);
%    mrkrData = zeros(nModelePoints,1);
%    xTData = zeros(nModelePoints,1); 
%% 降采样
    k=1;
    ki=1;
    frewind(fid);%将文件的位置指针移至开头部分
    kk=2;
    
    while ~feof(fid) 
        tline = fgetl(fid); % after reading each line, move onto the next line.
            if k > nPoints+1, break, end
            if k>1 & k==kk% starting from the second line.
                [tT xT yT zT] = strread(tline,'%f: %f %f %f',1);
                o_3D(folderIndex,ki,:) = [xT yT zT tT];
%                 o_3D(folderIndex,ki,1) = xT;
%                 o_3D(folderIndex,ki,2) = yT;
%                 o_3D(folderIndex,ki,3) = zT;
%                 o_3D(folderIndex,ki,4) = tT;
                ki=ki+1;
                kk=kk+incrt;
            end % of if k>1
            k=k+1;
    end % of while statement
%%
    initial_time=o_3D(folderIndex,1,4);
    for j=1:nPoints % make time start from zero
        o_3D(folderIndex,j,4)=o_3D(folderIndex,j,4)-initial_time;        
    end % of for      
    fclose(fid);
    % read Modeler.log file done...
end % for dbfolders

return