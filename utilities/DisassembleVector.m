function [ a1, a2 ] = DisassembleVector( vector, unknown_unfixed, known_unfixed)
%DISASSEMBLEMATRIX divides matrix in to parts corresponding to prescribed
% and unprescribed dofs and coupling parts:
% a1           corresponding to unprescribed, free dofs
% a2           corresponding to prescribed, free dofs

 a1=vector(unknown_unfixed);
 a2=vector(known_unfixed);

end

