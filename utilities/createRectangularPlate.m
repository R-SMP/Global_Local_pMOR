function [ plate, x0, xl, y0, yl ] = createRectangularPlate( lx, ly, nx, ny, varargin )
%CREATERECTANGULARPLATE Creates the mesh of a rectangular plate
%   inputs:
%       lx: length in x
%       ly: length in y
%       nx: number of elements in x
%       ny: number of elements in y
%       'elementType': used element type
%    outputs:
%        plate: FemModel with all nodes and elements
%        x0: nodes on x=0
%        xl: nodes on x=ly
%        y0: nodes on y=0
%        yl: nodes on y=lx

p = inputParser();
p.addParameter('elementType','ReissnerMindlinElement3d4n');
p.parse(varargin{:});
elementType = p.Results.elementType;

plate = FemModel();
nNodes = str2double(elementType(end-1));

%create nodes
idn = 1;
if any([3 4] == nNodes)
    for idy = 1:(ny+1)
        y = (idy-1) * ly / ny;
        for idx = 1:(nx+1)
            x = (idx-1) * lx / nx;
            plate.addNewNode(idn,x,y,0);
            idn = idn + 1;
        end
    end
elseif 8 == nNodes
    for idy = 1:(2*ny+1)
        y = (idy-1) * ly / (2*ny);
        if mod(idy,2)
            for idx = 1:(2*nx+1)
                x = (idx-1) * lx / (2*nx);
                plate.addNewNode(idn,x,y,0);
                idn = idn + 1;
            end
        else
            for idx = 1:(nx+1)
                x = (idx-1) * lx / nx;
                plate.addNewNode(idn,x,y,0);
                idn = idn + 1;
            end
        end
    end
elseif any([6 9] == nNodes)
    for idy = 1:(2*ny+1)
        y = (idy-1) * ly / (2*ny);
        for idx = 1:(2*nx+1)
            x = (idx-1) * lx / (2*nx);
            plate.addNewNode(idn,x,y,0);
            idn = idn + 1;
        end
    end
else
    msg = ['Utility createRectangularPlate: Elements with ', num2str(nNodes), ...
        ' nodes can not be used to mesh a rectangular plate.'];
    e = MException('MATLAB:bm_mfem:unsupportedElement',msg);
    throw(e);
end

%create elements
ide = 1;
if nNodes == 3
    for idy = 1:ny
        for idx = 1:nx
            plate.addNewElement(elementType,ide, ...
                [idx+(idy-1)*(nx+1) idx+(idy-1)*(nx+1)+1 idx+idy*(nx+1)]);
            plate.addNewElement(elementType,ide+1, ...
                [idx+(idy-1)*(nx+1)+1 idx+idy*(nx+1)+1 idx+idy*(nx+1)]);
            ide = ide + 2;
        end
    end
elseif nNodes == 4
    for idy = 1:ny
        for idx = 1:nx
            plate.addNewElement(elementType,ide, ...
                [idx+(idy-1)*(nx+1) idx+(idy-1)*(nx+1)+1 idx+idy*(nx+1)+1 idx+idy*(nx+1)]);
            ide = ide + 1;
        end
    end
elseif nNodes == 6
    idm = 0;
    for idy = 1:ny
        idu = idm(end)+1:idm(end)+2*nx+1;
        idm = idu(end)+1:idu(end)+2*nx+1;
        ido = idm(end)+1:idm(end)+2*nx+1;
        for idx = 1:nx
            plate.addNewElement(elementType,ide, ...
                [idu(idx*2-1) idu(idx*2+1) ido(idx*2-1) ...
                idu(idx*2) idm(idx*2) idm(idx*2-1)]);
            plate.addNewElement(elementType,ide+1, ...
                [idu(idx*2+1) ido(idx*2+1) ido(idx*2-1) ...
                idm(idx*2+1) ido(idx*2) idm(idx*2)]);
            ide = ide + 2;
        end
    end
elseif nNodes == 8
    idm = 0;
    x0id = [];
    xlid = [];
    for idy = 1:ny
        idu = idm(end)+1:idm(end)+2*nx+1;
        idm = idu(end)+1:idu(end)+nx+1;
        ido = idm(end)+1:idm(end)+2*nx+1;
        x0id = [x0id idu(1) idm(1) ido(1)]; %#ok<AGROW>
        xlid = [xlid idu(end) idm(end) ido(end)]; %#ok<AGROW>
        for idx = 1:nx
            plate.addNewElement(elementType,ide, ...
                [idu(idx*2-1) idu(idx*2+1) ido(idx*2+1) ido(idx*2-1) ...
                idu(idx*2) idm(idx+1) ido(idx*2) idm(idx)]);
            ide = ide + 1;
        end
    end
elseif nNodes == 9
    idm = 0;
    for idy = 1:ny
        idu = idm(end)+1:idm(end)+2*nx+1;
        idm = idu(end)+1:idu(end)+2*nx+1;
        ido = idm(end)+1:idm(end)+2*nx+1;
        for idx = 1:nx
            plate.addNewElement(elementType,ide, ...
                [idu(idx*2-1) idu(idx*2+1) ido(idx*2+1) ido(idx*2-1) ...
                idu(idx*2) idm(idx*2+1) ido(idx*2) idm(idx*2-1) idm(idx*2)]);
            ide = ide + 1;
        end
    end
end

%set boundaries
if any([3 4] == nNodes)
    nNodes = (nx+1)*(ny+1);
    x0 = plate.getNodes(1:nx+1:nNodes);
    y0 = plate.getNodes(1:nx+1);
    xl = plate.getNodes(nx+1:nx+1:nNodes);
    yl = plate.getNodes(nNodes-nx:nNodes);
elseif any([6 9] == nNodes)
    nNodes = (nx*2+1)*(ny*2+1);
    x0 = plate.getNodes(1:nx*2+1:nNodes);
    y0 = plate.getNodes(1:nx*2+1);
    xl = plate.getNodes(nx*2+1:nx*2+1:nNodes);
    yl = plate.getNodes(nNodes-nx*2:nNodes);
elseif nNodes == 8
    nNodes = ido(end);
    x0 = plate.getNodes(unique(x0id));
    y0 = plate.getNodes(1:2*nx+1);
    xl = plate.getNodes(unique(xlid));
    yl = plate.getNodes(nNodes-2*nx:nNodes);
end

end

