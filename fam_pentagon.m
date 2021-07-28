% function identifies pentagon-region, normalizes pentagon-coordinates

% input: alley & pentagon xy-coordinates
% output: pentagon xy-coordinates, pentagon-polyshapes

function [cP_x,cP_y,cP,pentagon_x,pentagon_y, inPentagon,...
    outPentagon] = fam_pentagon(alley_x,alley_y,pentagon_x, pentagon_y,xmin,xmax,ymin,ymax)

% coordinates pentagon-corner (normalized!)

[row,pentagonCornerNo] = size(pentagon_x);
for pentagon=1:pentagonCornerNo
    pentagon_x(row,pentagon)=datanorm(pentagon_x(row,pentagon),xmin,xmax);
end

for pentagon=1:pentagonCornerNo
    pentagon_y(row,pentagon)=datanorm(pentagon_y(row,pentagon),ymin,ymax);
end

innerAlleys_x= [alley_x(4,1); alley_x(3,1); alley_x(4,2);alley_x(3,2);alley_x(4,3);alley_x(3,3);...
    alley_x(4,4);alley_x(3,4);alley_x(4,5);alley_x(3,5)];
innerAlleys_y= [alley_y(4,1);alley_y(3,1);alley_y(4,2);alley_y(3,2);alley_y(4,3);alley_y(3,3);...
    alley_y(4,4);alley_y(3,4);alley_y(4,5);alley_y(3,5)];

inPentagon = [pentagon_x; pentagon_y]';
outPentagon = [innerAlleys_x, innerAlleys_y];

cP_x=[alley_x(4,1);alley_x(3,1);alley_x(4,2);alley_x(3,2);alley_x(4,3);alley_x(3,3);...
    alley_x(4,4);alley_x(3,4);alley_x(4,5);alley_x(3,5);NaN;...
    pentagon_x(1,1);pentagon_x(1,5);pentagon_x(1,4);pentagon_x(1,3);pentagon_x(1,2);pentagon_x(1,1)];
cP_y=[alley_y(4,1);alley_y(3,1);alley_y(4,2);alley_y(3,2);alley_y(4,3);alley_y(3,3);...
    alley_y(4,4);alley_y(3,4);alley_y(4,5);alley_y(3,5);NaN;...
    pentagon_y(1,1);pentagon_y(1,5);pentagon_y(1,4);pentagon_y(1,3);pentagon_y(1,2);pentagon_y(1,1)];
cP=polyshape(cP_x,cP_y);
end