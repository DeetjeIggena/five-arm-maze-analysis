% last @update 210527
% function assigns group-names according to ID
% adjust names as desired

function [group,Group]=fam_callGroup(ID)

% determine group, provide group names as string
if str2double(ID(2)) == 2 % group 0=Control, 1=Experimentalgroup
    group = 1;
    Group = 'g_1';
elseif str2double(ID(2)) == 4 % group 0=Control, 1=Experimentalgroup
    group = 2;
    Group = 'g_2';
elseif str2double(ID(2)) == 3
    group = 0;
    Group = 'Control';
else
    group=100;
    Group='F';
    disp('The provided id is out of limits. Group is set to "F"')
end
end
