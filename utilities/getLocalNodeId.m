function localNodeId = getLocalNodeId(nodeArray, x,y)
% getLocalNodeId: search for local node Id through given coordinates
% nodeArray: node array containing searched node
% x,y,z: coordinates of searched node

xCoor = nodeArray.getX;
yCoor = nodeArray.getY;


localNodeId_x = find((x-xCoor) == min(abs(x-xCoor)));

y_considered = zeros(length(localNodeId_x),1);

for j=1:length(localNodeId_x)
    y_considered(j) = yCoor(localNodeId_x(j));
end

id = find((y-y_considered) == min(abs(y-y_considered)));
localNodeId = localNodeId_x(id);

end


    
