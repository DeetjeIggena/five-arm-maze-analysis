clear; close all; clc; format compact;

%% five-arm-maze-Analysis
% @ date 200923 @ author Deetje Iggena
% @ date 210527 last update
% analysis-script for maze version 190819
% Matlab R2020b

% The script requires .csv as input-files.
% Input-file "trial_results" should contain info about the order of
% track-files and goal/-task information

% Input-files "trackepoint_movment_filename" should contain
% timestamp ("time"), x-position (pos_x), y-position (pos_z)

% BE AWARE:
% In case you change the input-format/ structure of input-files,
% please check which coloumns do
% contain your data and adjust the script accordingly

% Script starts with input:
%   1. Provide the span of subjects you would like to analyse, start with
%   the first subjectnumber, end with the last subjectnumber (Attention:
%   the relevant folder  (e.g. S001 or S002) has to be in a folder that is named after the
%   subjectnumber), the script will save the existing IDs into
%   group-arrays.
%   The group-arrays will be called later during the actual analysis

% optional: in case you would like to add participants later provide the
% length/ subjectnumber already analysed

%   2. results-folder & trackfolder will be created, choose your name

% Optional: if desired provide info about your projects to be saved in
% project-table

%   3. a m-table called ytable (your table) either exists (due to a previous analysis) 
%   or is created which will be saved in the folders provided
%   by you, the table will contain the relevant data

%   4. the actual analysis is embedded in three loops, 1st sessin loop (sN =
%   max number of sessions), 2nd group-loop (gN = max number of groups), 
%   3rd participant loop (called from group-arrays)

%   Optional: add analysis functions into loop

%   5. data is saved in m-table

% Quit by entering CTRL+C

%% Provide folder information
currentDirectory        = pwd; % contains data-folder
addpath(genpath(currentDirectory)); % add subfolders containing functions
finalFolderString       = 'S001'; % default
gN                      = 3; % default --> no of groups
sN                      = 2; % default
subjectsAnalysedExp     = 0;
subjectsAnalysedExp_2   = 0;
subjectsAnalysedCon     = 0;

%% select participants

[subject_start,subject_end]=sm_inputSubjectsNo();

% check for existing participants, alternativeley provide participant-list/array
% sort participants into group-arrays
e = 1; c = 1; d = 1;
subArrayExp = []; subArrayCon = []; subArrayExp_2 = [];

for sub = subject_start:subject_end
    subString = num2str(sub);
    folderIn = [currentDirectory '\' subString '\' finalFolderString];
    if exist(folderIn)
        g = sscanf(subString(2), '%d');
        if g == 2
            subArrayExp(e,1) = sub;
            e = e + 1;
        elseif g == 4
            subArrayExp_2(d,1) = sub;
            d = d + 1;
        else
            subArrayCon(c,1) = sub;
            c = c + 1;
        end
    end
end


%% create result & track-folder
% construct a path to results-table
resultFolder=[currentDirectory '\yourResultFolder']; %provide the name for your result folder
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
    disp('Your folder didn''t exist, a new result-folder was created')
end

%% create table for saving data
% construct the target-file/ load existing target-table
targetFileName         = '\yourTable.mat';
targetFilePath         = fullfile(resultFolder, targetFileName);
if isfile(targetFilePath)
    load(targetFilePath, 'ytable');
end
% initialize table of results if not-existing
if ~exist('table','var')
    ytable = [];
    save(targetFilePath, 'ytable');
end

%% Create starmaze
% Min-Max-values
values = table2array(readtable('ymaze_valuesMinMax.csv'));
ytable.xmin = values(1,1); ytable.xmax = values(2,1);
ytable.ymin = values(1,2); ytable.ymax = values(2,2);

% start-positions
start = table2array(readtable('ymaze_start.csv'));
[~,~,ytable.startPosition] = fam_normalizePosition(start,ytable.xmin,ytable.xmax,ytable.ymin,ytable.ymax);
% goal-positions
goal = table2array(readtable('ymaze_goal.csv'));
[ytable.goal_x,ytable.goal_y, ~] = fam_normalizePosition(goal,ytable.xmin,ytable.xmax,ytable.ymin,ytable.ymax);

% coordinates alley-corner
alley_x=table2array(readtable('ymaze_alley_x.csv'));
[cornerNo,alleyNo] = size(alley_x);
for alley=1:alleyNo
    for corner=1:cornerNo
        alley_x(corner,alley) = datanorm(alley_x(corner,alley),ytable.xmin,ytable.xmax);
    end
end
alley_y=table2array(readtable('ymaze_alley_y.csv'));
for alley=1:alleyNo
    for corner=1:cornerNo
        alley_y(corner,alley) = datanorm(alley_y(corner,alley),ytable.ymin,ytable.ymax);
    end
end

% combined pentagon
pentagon_x = table2array(readtable('ymaze_pentagon_x.csv'));
pentagon_y = table2array(readtable('ymaze_pentagon_y.csv'));
[cP_x,cP_y,cP,pentagon_x,pentagon_y, inPentagon, outPentagon] = fam_pentagon(alley_x,alley_y,pentagon_x, pentagon_y,ytable.xmin,...
    ytable.xmax,ytable.ymin,ytable.ymax);

% alley_polyshape
[alley_full_x,alley_full_y,ytable.alley_polyshape, alley_half_out_x, alley_half_out_y, alley_polyshape_1,...
    alley_half_in_x, alley_half_in_y, alley_polyshape_2] = fam_alleyPolyshape(alley_x,alley_y);
% rectangle_polyshape
[rec_x,rec_y,ytable.rec] = fam_rectPolyshape(alleyNo,alley_x, alley_y, pentagon_x, pentagon_y);
% triangle_polyshape
[tri_x,tri_y,ytable.tri] = fam_trianglePolyshape(alleyNo,alley_x, alley_y, pentagon_x, pentagon_y);
% array contains alley & central pentagon polyshapes
polyshape_array=[alley_polyshape_1{1,1} alley_polyshape_2{1,1} alley_polyshape_1{1,2} alley_polyshape_2{1,2}...
    alley_polyshape_1{1,3} alley_polyshape_2{1,3} alley_polyshape_1{1,4} alley_polyshape_2{1,4}...
    alley_polyshape_1{1,5} alley_polyshape_1{1,5} alley_polyshape_2{1,5} cP];

ytable.goalAlley = 0;
for i=1:length(alley_full_x)
    if inpolygon(ytable.goal_x, ytable.goal_y, alley_full_x(:,i), alley_full_y(:,i))
        ytable.goalAlley = i;
    end
end

% Test-Figure
% sm_wp3_testfigure(polyshape_array,ytable.goal_x,ytable.goal_y,...
%     ytable.start_x,ytable.start_y);

% -------------------------------------------------------------------------
%% Analysis
% -------------------------------------------------------------------------
for se = 1:sN % loop through sessions
    
    for gr = 1:gN % loop through groups
        if gr == 1
            g = 1;
            partArray = subArrayExp; % array experimental group
            s = subjectsAnalysedExp + 1;
        elseif gr == 2
            g = 2;
            partArray = subArrayExp_2; % array experimental group 2
            s = subjectsAnalysedExp_2 + 1;
        else
            g = gN;
            partArray = subArrayCon; % array control group
            s = subjectsAnalysedCon + 1;
        end
        
        participantLength   = length(partArray);
        
        for p = 1: participantLength
            
            partArray(p)
            ytable.session(se).group(g).subject(s).id             = partArray(p);
            ytable.session(se).group(g).subject(s).idString       = num2str(partArray(p));
            [ytable.session(se).group(g).subject(s).groupNo,...
                ytable.session(se).group(g).subject(s).groupName] = fam_callGroup(ytable.session(se).group(g).subject(s).idString); % group info
            
            finalFolderString     = [ 'S00' num2str(se)]; % default
            folderIn = [currentDirectory '\' ytable.session(se).group(g).subject(s).idString '\' finalFolderString];
            
            files = dir(fullfile(folderIn,'*.csv'));
            files = {files.name};
            
            logIndex = find(contains(files,'log'));
            files(:,logIndex)   = [];
            
            trialIndex = find(contains(files,'trial_results.csv'));
            files(:,trialIndex) = [];
            
            data_trial = readtable([folderIn, '\trial_results.csv']); % read in trial-file-info
            
            k     = 1; trial = 1;
            
            ytable.session(se).group(g).subject(s).guidance = 0;
            for f = 1:numel(files)
                
                guidance = strcmp(data_trial.trial_type(f,1),'guidance');
                
                if guidance == 1
                    ytable.session(se).group(g).subject(s).guidance = 1; % exclude guidance-trial from analysis
                    continue
                end
                
                name = files{f};
                name = fullfile(folderIn, name);
                data = readtable(name);
                
                %% data preparation -> remove pause, first two rows, empty data
                % --------------------------------------------------------------------------------------
                if isempty(data) % check whether file contains any data
                    trial = trial + 1;
                    continue
                end
                
                data(1:2,:)     = [];
                j = 1; data_paused = [];
                for i=1:length(data.time)-1
                    samplingRate(i)    = data.time(i+1,1) - data.time(i,1);
                    if strcmp(data.gameIsPaused(i),'True')
                        data_paused(j,1) = data.time(i);
                        j = j+1;
                    end
                end

                avgSamplingRate        = sum(samplingRate)/length(samplingRate);

                data(strcmp(data.gameIsPaused,'True'),:) = []; % remove xy during pause
                if isempty(data)
                    trial = trial + 1;
                    continue
                end
                
         
                %% Info depending on single trial: Feedback, Trial-Type/Condition, goal-& startposition
                % --------------------------------------------------------------------------------------
                
                ytable.session(se).group(g).subject(s).sm.trial(k).trial           = trial;
                ytable.session(se).group(g).subject(s).sm.trial(k).block_num       = data_trial.block_num(f,1);
                ytable.session(se).group(g).subject(s).sm.trial(k).trial_num       = data_trial.trial_num(f,1);
                ytable.session(se).group(g).subject(s).sm.trial(k).trial_in_block  = data_trial.trial_num_in_block(f,1);
                ytable.session(se).group(g).subject(s).sm.trial(k).trial_type      = data_trial.trial_type(f,1);
                ytable.session(se).group(g).subject(s).sm.trial(k).feedback        = data_trial.trial_feedback(f,1);

                    
                %% select data
                t = data.time; x = data.pos_x; y = data.pos_z;
                
                x = datanorm(x,ytable.xmin,ytable.xmax);
                y = datanorm(y,ytable.ymin,ytable.ymax);   % data normalization for coordinates
                
                data_length = length(x);
                
                % save coordinates in table for further pr
                ytable.session(se).group(g).subject(s).sm.trial(k).x         = x;
                ytable.session(se).group(g).subject(s).sm.trial(k).y         = y;
                
                %%  Variables depending on starting-positions
                % --------------------------------------------------------------------------------------
                % determine start
                [ytable.session(se).group(g).subject(s).sm.trial(k).start]   = fam_trialStart(ytable.session(se).group(g).subject(s).sm.trial(k).trial_type);
                
                [nodesNo, nodesArray] = fam_shortestPath(ytable.startPosition,...
                    inPentagon, outPentagon, ytable.session(se).group(g).subject(s).sm.trial(k).start, ytable.goalAlley);
                
                %% data-analysis
                % --------------------------------------------------------------------------------------
                % Time-Analysis using timestamp
                b = t(end,1); a = t(1,1);
                
                [ytable.session(se).group(g).subject(s).sm.trial(k).time,...
                    ~] = fam_time(a,b,data_length,avgSamplingRate);
                
                % Coordinate-Analysis using x & y, z
                % -------------------------------------------------------------
                % Path % distance-Analysis
                pathLength       = zeros(1, data_length); dist_to_goal     = zeros(1, data_length);
                
                for i=1:data_length-1
                    pathLength(i)   = fm_distance(x(i),x(i+1),y(i),y(i+1));% cumulative distance traveled
                    dist_to_goal(i) = fam_distance(x(i),ytable.goal_x,y(i),ytable.goal_y); % cumulative distance to allocentric target
                end
                
                ytable.session(se).group(g).subject(s).sm.trial(k).pathLength         = sum(pathLength); % path-length/ distance travelled
                ytable.session(se).group(g).subject(s).sm.trial(k).avgDistTarget      = sum(dist_to_goal)/data_length; % cumulative distance to allocentric target
                ytable.session(se).group(g).subject(s).sm.trial(k).finalDistance      = fam_distance(x(end,:),ytable.goal_x,y(end,:),ytable.goal_y); % final distance to target

                % Velocity
                ytable.session(se).group(g).subject(s).sm.trial(k).velocity       = ytable.session(se).group(g).subject(s).sm.trial(k).pathLength/ytable.session(se).group(g).subject(s).sm.trial(k).time;
                                               
                %% Zone-Analysis
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry] = fam_coordinatesArea(x,y,alley_full_x,alley_full_y,ytable.session(se).group(g).subject(s).sm.trial(k).time);
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_abs_all    = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_abs);
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel_all    = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel);
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time_all   = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time);
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry_all  = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry);
                
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_rel, ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_time,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_entry] = fam_coordinatesArea(x,y,rec_x,rec_y, ytable.session(se).group(g).subject(s).sm.trial(k).time);
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_triangle_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_triangle_rel, ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_triangle,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_triangle_entry]  = fam_coordinatesArea(x,y,tri_x,tri_y, ytable.session(se).group(g).subject(s).sm.trial(k).time);
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_pentagon_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_pentagon_rel, ytable.session(se).group(g).subject(s).sm.trial(k).zone_pentagon_time,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_pentagon_entry]  = fam_coordinatesArea(x,y,tri_x,tri_y, ytable.session(se).group(g).subject(s).sm.trial(k).time);
                
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_rel, ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_time,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_entry]     = fam_coordinatesArea(x,y,alley_half_out_x,alley_half_out_y, ytable.session(se).group(g).subject(s).sm.trial(k).time);
                [ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_abs,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_rel, ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_time,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_entry]      = fam_coordinatesArea(x,y,alley_half_in_x,alley_half_in_y, ytable.session(se).group(g).subject(s).sm.trial(k).time);
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_entries_all      = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_entry) + ...
                    sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_entry) + sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_entry) + ...
                    sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_triangle_entry);
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zone_entries_allRec   = sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyOut_entry) + ...
                    sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_alleyIn_entry) + ...
                    sum(ytable.session(se).group(g).subject(s).sm.trial(k).zone_rectangle_entry);
                
                if ytable.goalAlley == ytable.session(se).group(g).subject(s).sm.trial(k).start
                    goalAlleyEntry      = 0;
                    goalAlleyTime       = 0;
                    goalAlleyRel        = 0;
                else
                    goalAlleyEntry      = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry(ytable.goalAlley);
                    goalAlleyTime       = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time(ytable.goalAlley);
                    goalAlleyRel        = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel(ytable.goalAlley);
                end
                
                if ytable.session(se).group(g).subject(s).sm.trial(k).start < 6
                    startAlleyEntry     = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry(ytable.session(se).group(g).subject(s).sm.trial(k).start);
                    startAlleyTime      = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time(ytable.session(se).group(g).subject(s).sm.trial(k).start);
                    startAlleyRel       = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel(ytable.session(se).group(g).subject(s).sm.trial(k).start);
                else
                    startAlleyEntry      = 0;
                    startAlleyTime       = 0;
                    startAlleyRel        = 0;
                end
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyEntryCorrect   = startAlleyEntry + ...
                    goalAlleyEntry;
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyTimeCorrect    = startAlleyTime + ...
                    goalAlleyTime;
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyRelCorrect     = startAlleyRel + ...
                    goalAlleyRel;
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyRelIncorrect   = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_rel_all - ...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyRelCorrect;
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyTimeIncorrect  = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_time_all - ...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyTimeCorrect;
                
                ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyEntryIncorrect = ytable.session(se).group(g).subject(s).sm.trial(k).zone_alley_entry_all - ...
                    ytable.session(se).group(g).subject(s).sm.trial(k).zoneAlleyEntryCorrect;
                
                %% Exploration-Analysis
                % determine final position --> alley-number or pentagon
                [ytable.session(se).group(g).subject(s).sm.trial(k).finalAlley,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).finalAlleyNo,...
                    ytable.session(se).group(g).subject(s).sm.trial(k).final_pentagon] = fam_finalZone(x,y,alley_full_x,alley_full_y,cP_x, cP_y);
                
                % determine whether final position is in the correct final alley
                if ytable.session(se).group(g).subject(s).sm.trial(k).finalAlleyNo == ytable.goalAlley
                    ytable.session(se).group(g).subject(s).sm.trial(k).correctAlley = 1;
                else
                    ytable.session(se).group(g).subject(s).sm.trial(k).correctAlley = 0;
                end
                
                % success
                dis_criterium = 0.1; % default - criterium
                ytable.session(se).group(g).subject(s).sm.trial(k).success     = fam_success(dis_criterium, ytable.session(se).group(g).subject(s).sm.trial(k).finalDistance);

                
                %% ------- Ideal measures & Error calculation -------------------------------------
                % calculation path-error/deviation/accuracy
                % calculation distance-error/deviation/accuracy
                
                % interpolate data for further analysis
                [nodesNum, ~] = size(nodesArray);

                if nodesNum > 1
                    [xi_al, yi_al]     = fam_interpolation(nodesArray,avgSamplingRate);
                    interp_length      = length(xi_al);
                    idealPathLength    = zeros(1, interp_length);
                    idealDist_to_goal  = zeros(1, interp_length);
                    
                    for i=1:interp_length-1
                        idealPathLength(i)    = fam_distance(xi_al(i),xi_al(i+1),yi_al(i),yi_al(i+1));% cumulative distance traveled
                        idealDist_to_goal(i)  = fam_distance(xi_al(i),ytable.goal_x,yi_al(i),ytable.goal_y); % cumulative distance to allocentric target
                    end
                    
                    ytable.session(se).group(g).subject(s).sm.trial(k).idealPathLength    = sum(idealPathLength); % path-length/ distance travelled
                    ytable.session(se).group(g).subject(s).sm.trial(k).pathError          = fam_error(ytable.session(se).group(g).subject(s).sm.trial(k).pathLength,...
                            ytable.session(se).group(g).subject(s).sm.trial(k).idealPathLength);
                    ytable.session(se).group(g).subject(s).sm.trial(k).idealAvgDistTarget = sum(idealDist_to_goal)/interp_length; % cumulative distance to allocentric target
                    ytable.session(se).group(g).subject(s).sm.trial(k).distanceError      = fam_error(ytable.session(se).group(g).subject(s).sm.trial(k).avgDistTarget,...
                            ytable.session(se).group(g).subject(s).sm.trial(k).idealAvgDistTarget);
                else

                    if ytable.session(se).group(g).subject(s).sm.trial(k).pathLength         <= 0.1
                        ytable.session(se).group(g).subject(s).sm.trial(k).idealPathLength    = 0;
                        ytable.session(se).group(g).subject(s).sm.trial(k).idealAvgDistTarget = 0;
                        ytable.session(se).group(g).subject(s).sm.trial(k).pathError          = 0;
                        ytable.session(se).group(g).subject(s).sm.trial(k).distanceError      = 0;
                    else 
                        ytable.session(se).group(g).subject(s).sm.trial(k).idealPathLength    = 0.05;
                        ytable.session(se).group(g).subject(s).sm.trial(k).pathError          = fam_error(ytable.session(se).group(g).subject(s).sm.trial(k).pathLength,...
                                ytable.session(se).group(g).subject(s).sm.trial(k).idealPathLength);
                        ytable.session(se).group(g).subject(s).sm.trial(k).idealAvgDistTarget = 0.025;
                        ytable.session(se).group(g).subject(s).sm.trial(k).distanceError      = fam_error(ytable.session(se).group(g).subject(s).sm.trial(k).avgDistTarget,...
                                ytable.session(se).group(g).subject(s).sm.trial(k).idealAvgDistTarget);
                    end
                end

                %% Save data
                % ---------------------------------------------------------------------------------------------------------------------------
                % construct the target-file/ load existing target-table
                save(targetFilePath, 'ytable','-append');
                
                clear data x y r t xi_al yi_al goalPos data_paused samplingRate
                %%
                k     = k +1;
                trial = trial + 1;
            end
            clear data_trial logIndex trialIndex folderIn name
            s = s + 1;
        end
    end
end

% %---------------------------------------------------------------------------------------------------------------------------
% %% Save data
% % ---------------------------------------------------------------------------------------------------------------------------
% % construct the target-file/ load existing target-table
% save(targetFilePath, 'ytable', '-append');
% 
% close all
% clear variables
