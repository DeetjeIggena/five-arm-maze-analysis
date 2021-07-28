% function creates nodes which consists of corners within the maze,
% according to the connectins and weights the shortest path between nodes
% (goal & start)is identified

% input = nodeCoordinates
% output = number of nodes and array of nodes included in the shortest path
% to target

% required to calculate ideal measures

function [nodesNo, nodesArray] = fam_shortestPath(arrayStart,...
    arrayInPentagon, arrayExPentagon, startPos, finalPos)
if startPos == 2
    startPos = 4;
end
%%
% 20 nodes by default
smNodes = [1 1 1 2 2 2 3 3 3 ...
    4 4 4 5 5 5 ... 
    6 6 6 6 7 7 7 7 8 8 8 8 ...
    14 16 18 20 22 ...
    9 10 11 12 13 ...
    10 11 11 12];
smConnection = [23 9 14 15 10 16 17 11 18 ...
    19 12 20 21 13 22 ...
    14 15 9 10 12 13 20 21 9 13 22 23 ...
    15 17 19 21 23 ...
    10 11 12 13 9 ...
    17 16 19 18];

smWeights = [3 4 3 3 4 3 3 4 3 3 4 3 3 4 3 ...
    1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 ...
    2 2 2 2 2 ...
    2 2 2 2 2 ...
    2 2 2 2];

x = [arrayStart(1,1) arrayStart(2,1) arrayStart(3,1) arrayStart(4,1) arrayStart(5,1) ...
    arrayStart(6,1) arrayStart(7,1) arrayStart(8,1) ...
    arrayInPentagon(1,1) arrayInPentagon(2,1) arrayInPentagon(3,1)...
    arrayInPentagon(4,1) arrayInPentagon(5,1) ...
    arrayExPentagon(2,1) arrayExPentagon(3,1) arrayExPentagon(4,1) ...
    arrayExPentagon(5,1) arrayExPentagon(6,1) arrayExPentagon(7,1) ...
    arrayExPentagon(8,1) arrayExPentagon(9,1) arrayExPentagon(10,1) ...
    arrayExPentagon(1,1) ];
y = [arrayStart(1,2) arrayStart(2,2) arrayStart(3,2) arrayStart(4,2) arrayStart(5,2) ...
    arrayStart(6,2) arrayStart(7,2) arrayStart(8,2)...
    arrayInPentagon(1,2) arrayInPentagon(2,2) arrayInPentagon(3,2)...
    arrayInPentagon(4,2) arrayInPentagon(5,2) ...
    arrayExPentagon(2,2) arrayExPentagon(3,2) arrayExPentagon(4,2) ...
    arrayExPentagon(5,2) arrayExPentagon(6,2) arrayExPentagon(7,2) ...
    arrayExPentagon(8,2) arrayExPentagon(9,2) arrayExPentagon(10,2) ...
    arrayExPentagon(1,2)];

nodesCoordinates = [x;y]';

G = graph(smNodes,smConnection,smWeights);
names = {'S1' 'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8' ...
    'iP1' 'iP2' 'iP3' 'iP4' 'iP5' ...
    'eP1' 'eP2' 'eP3' 'eP4' 'eP5' 'eP6' 'eP7' 'eP8' 'eP9' 'eP10'}';
G.Nodes.Name = names;

nodesNo = length(names);

figure('visible','off')
 p = plot(G,'XData',x,'YData',y);

[pathNodes,~] = shortestpath(G,startPos,finalPos);
highlight(p,pathNodes,'EdgeColor','g');

[~,col] = size(pathNodes);

%nodesCoordinates
for i=1:col
    z = pathNodes(i);
    nodesArray(i,1) = nodesCoordinates(z,1);
    nodesArray(i,2) = nodesCoordinates(z,2);
    
end

end