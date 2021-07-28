% function returns the starting-position
% adjust according to your paradigm

function [start]=fam_trialStart(trial_type)


start_1  = contains(trial_type,'training');
start_11 = contains(trial_type,'ego');
start_2  = contains(trial_type,'startC');
start_3  = contains(trial_type,'startE');
start_4  = contains(trial_type,'startG');
start_5  = contains(trial_type,'startI');
start_6  = contains(trial_type,'startAC');
start_7  = contains(trial_type,'startGI');
start_8  = contains(trial_type,'startIA');

if start_1 == 1 || start_11 == 1
    start = 1;
elseif start_6 == 1
    start = 6;
elseif start_7 == 1
    start = 7;
elseif start_8 == 1
    start = 8;
elseif start_2 == 1
    start = 2;
elseif start_3 == 1
    start = 3;
elseif start_4 == 1
    start = 4;
elseif start_5 == 1
    start = 5;
end

end