% function interpolates points between nodes, distance according to average
% sampling rate
% 'makima' was chosen to account for a biological person navigating

function [x_line, y_line] = fam_interpolation(nodesArray,samplingRate)

nodeFinal = length(nodesArray);
% interpolate data depending on start-positions 
if nodesArray(1) < nodesArray(nodeFinal)
    x_line = nodesArray(1):samplingRate:nodesArray(nodeFinal);
else
   x_line = nodesArray(nodeFinal):samplingRate:nodesArray(1);
end 
    
x_line = transpose(x_line);
y_line = interp1(nodesArray(:,1),nodesArray(:,2),x_line,'makima');

end
