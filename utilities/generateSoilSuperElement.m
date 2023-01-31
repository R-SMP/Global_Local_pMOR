function nodearray = generateSoilSuperElement(dx,Nx)
% create nodeArray of SoilSuperElement
x=0:dx:(Nx-1)*dx;
y=x;
z=0; j=0;

for ny = 1:length(y)
    for nx = 1:length(x)
        j=j+1;
        nodearray(j)=Node(j,x(nx),y(ny),z);
    end    
end

end