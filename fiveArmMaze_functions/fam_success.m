% function returns whether coordinates are in a certain distance to the
% target-position

% Input: distance to allo/ egocentric target-position
% Output: 0 = unsuccessful, 1 = successful

function [success]=fam_success(distance,final_distance)


if final_distance <= distance
    success=1;
else
    success=0;
end
 
end
