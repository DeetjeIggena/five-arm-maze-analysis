% function identifies triangle-region between central and peripheral
% alleys, required for detailed area-analysis

% input: number of alleys, alley & pentagon xy-coordinates
% output: triangle xy-coordinates, triangle-polyshapes

function [tri_x,tri_y,tri]=fam_trianglePolyshape(alleyNo,alley_x, alley_y, pentagon_x, pentagon_y)

% Define triangles
for alley=1:alleyNo
    tri_x(:,alley) = [alley_x(4,alley);alley_x(3,alley);pentagon_x(1,alley);alley_x(4,alley)];
    tri_y(:,alley) = [alley_y(4,alley);alley_y(3,alley);pentagon_y(1,alley);alley_y(4,alley)];
    tri{alley}     = polyshape(tri_x(:,alley),tri_y(:,alley));
end
end
