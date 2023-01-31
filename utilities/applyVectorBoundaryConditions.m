function [ vector ] = applyVectorBoundaryConditions( vector, bc )
%APPLYVECTORBOUNDARYCONDITIONS Removes bc entries from a vector
%   Detailed explanation goes here
if ~bc; return; end
    vector(bc) = [];

end

