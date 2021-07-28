% function identifies rectangle-region between central
% alleys & triangles, required for detailed area-analysis

% input: number of alleys, alley & pentagon xy-coordinates
% output: rectangle xy-coordinates, rectangle-polyshapes

function [rec_x,rec_y,rec] = fam_rectPolyshape(alleyNo,alley_x, alley_y, pentagon_x, pentagon_y)

for alley=1:alleyNo
    if alley==alleyNo
        x = -alleyNo;
    else
        x = 0;
    end
    
    rec_x(:,alley)  = [alley_x(3,alley);alley_x(4,alley+1+x);pentagon_x(1,alley+1+x);pentagon_x(1,alley);alley_x(3,alley)];
    rec_y(:,alley)  = [alley_y(3,alley);alley_y(4,alley+1+x);pentagon_y(1,alley+1+x);pentagon_y(1,alley);alley_y(3,alley)];
    rec{alley}      = polyshape(rec_x(:,alley),rec_y(:,alley));
    
end
end