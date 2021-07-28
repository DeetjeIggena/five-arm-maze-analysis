% function returns a testfigure for visual inspection of alley-arrangement
% and goal-position
% @update 210527

% input: array containing alley-ployshapes, goal-position
% output: figure

function fam_testfigure(polyshape_array,goal_x,goal_y, start_x, start_y)

figure('Position',[500 200 580 500]);
set(gca,'xtick',[0 1],'ytick',[0 1]);
plot(polyshape_array);
axis([0 1 0 1]);
title('Star-Maze');
hold on
for g=1:length(goal_x)
    viscircles([goal_x(g) goal_y(g)], 0.015);
end
for g=1:length(start_x)
    viscircles([start_x(g) start_y(g)], 0.005, 'Color','b');
end
hold off

text(0.1,0.8,'Alley 5');
text(0.9,0.9,'Alley 2');
text(0.1,0.2,'Alley 4');
text(0.9,0.2,'Alley 3');
text(0.6,0.9,'Alley 1');
end
