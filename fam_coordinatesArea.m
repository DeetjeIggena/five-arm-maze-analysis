% function returns number of coordinates, relative amount of coordinates in
% alley & entries into alley
% last @update 02.06.2021
% @autho D.Iggena

% Input: x- & y coordinates, area of interest, trial duration
% Output: number of coordinates, rel. amount of coordinates in area,
% entries into area

function [abs_zone,rel_zone, time_zone, entry]= fam_coordinatesArea(x,y,area_x,area_y, time)

[~,col]=size(area_x);
abs_zone=zeros(1,col); rel_zone=zeros(1,col); time_zone=zeros(1,col); entry=zeros(1,col);
data_length= length(x);

for c=1:col
    [in,~] = inpolygon(x,y,area_x(:,c),area_y(:,c));
    coordinates_in_area =numel(x(in));
    abs_zone(1,c) = abs_zone(1,c)+coordinates_in_area; %absolut
    rel_zone(1,c) = abs_zone(1,c)/data_length; % relativ
    time_zone(1,c) = time*rel_zone(1,c);
    
    % entries
    if inpolygon(x(1,1),y(1,1),area_x(:,c),area_y(:,c))
        entry(1,c)=entry(1,c)+1;
    end
    for i=2:data_length-1
        if inpolygon(x(i),y(i),area_x(:,c),area_y(:,c)) && ~inpolygon(x(i-1),y(i-1),area_x(:,c),area_y(:,c))
            entry(1,c)=entry(1,c)+1;
        end
    end
end
end
