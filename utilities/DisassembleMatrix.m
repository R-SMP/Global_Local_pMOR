function [ a11, a12, a21, a22 ] = DisassembleMatrix( matrix, unknown_unfixed, known_unfixed)
%DISASSEMBLEMATRIX divides matrix in to parts corresponding to prescribed
% and unprescribed dofs and coupling parts:
% a11           corresponding to unprescribed, free dofs
% a12, a21      corresponding to coupling parts
% a22           corresponding to prescribed, free dofs

 temp1=matrix(unknown_unfixed,:);
 temp2=matrix(known_unfixed,:);
 
 a11=temp1(:,unknown_unfixed);
 a12=temp1(:,known_unfixed);
 a21=temp2(:,unknown_unfixed);
 a22=temp2(:,known_unfixed);

end

