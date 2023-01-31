function [femModel] = appendUnitCells(femModel, n, eigenvalue, eigenvector)
%APPENDUNITCELLS appends n unit cells in x direction
%   Cell has to be periodic with respect to x!
%   EIGENVALUE and EIGENVECTOR can be specified in order to visualize wave
%   propagation using the Bloch theorem

%Get nodes and elements of original model
nodeArray = femModel.getAllNodes();
elementArray = femModel.getAllElements();

nNodes = length(nodeArray);
nElements = length(elementArray);

%Determine length of unit cell
l_x = min(nodeArray.getX());
r_x = max(nodeArray.getX());
lengthUnitCell = abs(r_x-l_x);

%Extend nodeArray
for i=1:n
    for j=1:nNodes
        x = nodeArray(j).getX;
        x = x + lengthUnitCell*i;
        y = nodeArray(j).getY;
        node = femModel.addNewNode(nNodes*i+j, x, y, 0);
        node.addDof(nodeArray(j).getDofMap().keys);
        fixed = nodeArray(j).getDofArray().isFixed();
        dofs = node.getDofArray();
        dofs(fixed).fix;
    end
end

%Extend elementArray
for i=1:n
    for j=1:nElements
        nodeIds = elementArray(j).getNodes().getId;
        nodeIds = nodeIds+i*nNodes;
        type = class(elementArray(j));
        femModel.addNewElement(type, nElements*i+j, nodeIds);
    end
end

femModel.updateDofArray();
femModel.setIsInitialized(true);


if nargin > 2
    for k = 1:length(eigenvalue)
        tmp = cell(n+1,1);
        tmp{1} = eigenvector(:,k);
        
        for i=2:n+1
            tmp{i} = eigenvalue(k)*tmp{i-1};
        end
        res = cell2mat(tmp);
        
        if k==1
            SimpleAssembler.assignResultsToDofs(femModel, res);
            femModel.getProcessInfo.setValue('FREQUENCY', eigenvalue(k));
        else
            SimpleAssembler.appendValuesToDofs(femModel, res);
            femModel.getProcessInfo.appendValue('FREQUENCY', eigenvalue(k));
        end
    end
end

end

