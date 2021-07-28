% Aim: compare the shape of two aligned trajectories. 

% input: x- & y-coordinates of two trajectories
% output: Dissimilarity-ratio, testfigure

% 1st block get data 
    % provide information about data location
    % skip step if your data is already organized in two matrices (2*n)

% 2nd block data analysis
    % script aligns two trajectories with different lenghts by dynamical time
    % warping (function dtw @mathworks)
    % script compares the shape of two trajectories with the same lenght by
    % procrustes analysis (function procrustes @mathworks)

% 3rd block plot figure

% @last update 210728 by D.Iggena

%% Get data
% Provide folder information
currentDirectory        = pwd; % contains data-folder
addpath(genpath(currentDirectory)); % add subfolders 
% load data
% construct the target-file/ load existing target-table
targetFileName    = '\your_file';
resultFolder      = [currentDirectory '\yourFolder']; 
targetFilePath    = fullfile(resultFolder, targetFileName);

load(targetFilePath, 'yourTable'); % table contains x & y coordinates, or provide coordinates from files such as .xlsx, .csv

% select trajectories for matrix
M_1 = [yourTable.x_1 yourTable.y_1]; M_2 = [yourTable.x_2 yourTable.y_2];
M1 = M_1'; M2 = M_2'; 
%% data-analysis

% align trajectories with dynamical time warping
[dtw_trajectory,i1,i2] = dtw(M1,M2);

% get selected coordinates from original trajectories
X  = [x_1(i1),y_1(i1)]; Y  = [x_2(i2),y_2(i2)];
            
% Procrustes
[dissimimlarityRatio,Z,~] = procrustes(X,Y,'Scaling',false, 'reflection',false); % set scaling and reflection according to your needs

%% Figure
figure
plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'y.',Z(:,1),Z(:,2),'b.');
