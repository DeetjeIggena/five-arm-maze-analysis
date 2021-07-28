% function to normalize data-point
% data normalization for coordinates

% input: datapoint, min & max
% output: normalized datapoint

function DN=datanorm(c,cmin,cmax)
DN=((c-cmin)/(cmax-cmin));

