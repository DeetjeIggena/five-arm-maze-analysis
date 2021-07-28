% function returns normalized position as matrix & vectors
% respectiveley
function [pos_x,pos_y,positions] = fam_normalizePosition(pos,xmin,xmax,ymin,ymax)

[row,~]=size(pos);

for r=1:row
    pos_x(r,1) = datanorm(pos(r,1),xmin,xmax);
    pos_y(r,1) = datanorm(pos(r,2),ymin,ymax);
end

positions = [pos_x, pos_y];

end
