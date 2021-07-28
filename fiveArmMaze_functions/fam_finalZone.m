% function returns final area/ alley or final position in pentagon

function [final_alley, final_alley_no, final_pentagon]=fam_finalZone(x,y,area_x,area_y,cP_x, cP_y)

[~,col]=size(area_x);
final_alley=zeros(1,col);

final_alley_no = 0;
for c=1:col
    final_alley(c) = inpolygon(x(end,:),y(end,:),area_x(:,c), area_y(:,c));
    if inpolygon(x(end,:),y(end,:),area_x(:,c), area_y(:,c))
        final_alley_no = c ;
    end
    final_pentagon = inpolygon(x(end,:),y(end,:),cP_x, cP_y);
    
end
