classdef GiDOutput < ModelIO
    %GIDOUTPUT Model output to GiD
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        
        function obj = GiDOutput(file)
            settings = struct('writeModelParts', false, ...
                'nodalResults', [], ...
                'elementalResults', []);
            
            if nargin == 0
                super_args = {};
            elseif nargin == 1
                super_args = {file, settings};
            end
            
            obj@ModelIO(super_args{:});
            
        end
        
        function p = getPreferences(obj)
            p = obj.settings;
        end
        
        function writeMesh(obj, model)
            fid = fopen([obj.file, '.post.msh'],'w');
            
            %write complete model
            elms = model.getAllElements;
            eTypes = string(arrayfun(@(x) class(x), elms, 'UniformOutput', 0));
            ueTypes = string(unique(eTypes));
            
            for ii=1:length(unique(eTypes))
                ce = elms(eTypes == ueTypes(ii));
                cName = class(ce);
                
                if isa(ce, 'DummyElement')
                    %dummy elements have to be treated differently
                    if length(ce) > 1; error('not supported yet'); end
                    elementTable = ce.getElementTable();
                    elementTable = string(elementTable(:,1));
                    [nc, et] = ce.getNodeConnectivity();
                    
                    for jj = 1:length(elementTable)
                        cAnsysElement = char(elementTable(jj));
                        [~, nnodes, etype] = AnsysInput.getAnsysElementInfo(cAnsysElement);
                        header = ['MESH "', cAnsysElement, '_Mesh_', ...
                            num2str(ii), '" dimension 3 ElemType ', etype, ...
                            ' Nnode ', num2str(nnodes), '\n'];
                        fprintf(fid, header);
                        
                        fprintf(fid, 'Coordinates\n');
                        cAnsysElementIds = et(et(:,2)==jj);
                        elementalNodeIds = cell2mat(nc(cAnsysElementIds))';
                        uniqueNodeIds = unique(reshape(elementalNodeIds,[],1));
                        uniqueNodes = model.getNodes(uniqueNodeIds);
                        coords = uniqueNodes.getCoords();
                        dat = sortrows([uniqueNodeIds coords], 1);
                        fprintf(fid, '%u %f %f %f \n', dat');
                        fprintf(fid, 'End Coordinates\n');
                        
                        if strcmp(cAnsysElement,'SOLID186')
                            %sort ANSYS node numbers according to GiD
                            elementalNodeIds = elementalNodeIds(:,...
                                [1 2 3 4 5 6 7 8 9 10 11 12 17 18 19 20 13 14 15 16]);
                        end
                        
                        fprintf(fid, 'Elements\n');
                        dat = [cAnsysElementIds elementalNodeIds];
                        format = ['%u ' repmat('%u ', 1, size(dat,2)-1) '\n'];
                        fprintf(fid, format, dat');
                        fprintf(fid, 'End Elements\n');
                    end
                    
                elseif isa(ce, 'SoilSuperElement')
                    %output for the SoilSuperElement
                    if length(ce) > 1; error('not supported yet'); end
                    
                    header = ['MESH "', cName, '_Mesh_', num2str(ii), ...
                        '" dimension 3 ElemType Quadrilateral Nnode 4\n'];
                    fprintf(fid, header);
                    
                    fprintf(fid, 'Coordinates\n');
                    nodes = unique(ce.getNodes());
                    ids = nodes.getId();
                    coords = nodes.getCoords();
                    dat = sortrows([ids' coords], 1);
                    fprintf(fid, '%u %f %f %f \n', dat');
                    fprintf(fid, 'End Coordinates\n');
                    
                    fprintf(fid, 'Elements\n');
                    nc = ce.getNodeConnectivity(); %nc is an n x 5 array with [sub element id node_id_1 2 3 4] 
                    fprintf(fid, '%u %u %u %u %u \n', nc');
                    fprintf(fid, 'End Elements\n');
                    
                else
                    %write standard elements
                    header = ['MESH "', cName, '_Mesh_', num2str(ii), ...
                        '" dimension ', num2str(cName(end-3)), ' ElemType ', ...
                        ce.getElementType, ' Nnode ', num2str(cName(end-1)), '\n'];
                    fprintf(fid, header);

                    fprintf(fid, 'Coordinates\n');
                    nodes = ce.getNodes();
                    uniqueNodes = unique(nodes);
                    ids = uniqueNodes.getId();
                    coords = uniqueNodes.getCoords();
                    dat = sortrows([ids' coords], 1);
                    fprintf(fid, '%u %f %f %f \n', dat');
                    fprintf(fid, 'End Coordinates\n');

                    fprintf(fid, 'Elements\n');
                    ids = ce.getId();
                    nodeIds = arrayfun(@(x) x.getId(), nodes);
                    dat = sortrows([ids' nodeIds], 1);
                    format = ['%u ' repmat('%u ', 1, size(dat,2)-1) '\n'];
                    fprintf(fid, format, dat');
                    fprintf(fid, 'End Elements\n');                    
                end
            end
            
            %write model parts
            if obj.settings.writeModelParts
                nelms = length(elms);
                mps = model.getAllModelParts();
                for ii=1:length(mps)
                    mp = model.getModelPart(mps{ii});
                    
                    %write model part elements
                    elms = mp.getElements();
                    eTypes = string(arrayfun(@(x) class(x), elms, 'UniformOutput', 0));
                    ueTypes = string(unique(eTypes));
                    
                    for jj=1:length(unique(eTypes))
                        ce = elms(eTypes == ueTypes(jj));
                        if isa(ce, 'DummyElement') || isa(ce, 'SoilSuperElement'); continue; end
                        cName = class(ce);
                        header = ['MESH "', mps{ii}, '_', cName, '_Mesh_', num2str(jj), ...
                            '" dimension ', num2str(cName(end-3)), ' ElemType ', ...
                            ce.getElementType, ' Nnode ', num2str(cName(end-1)), '\n'];
                        fprintf(fid, header);
                        
                        fprintf(fid, 'Coordinates\n');
                        nodes = ce.getNodes();
                        uniqueNodes = unique(nodes);
                        ids = uniqueNodes.getId();
                        coords = uniqueNodes.getCoords();
                        dat = sortrows([ids' coords], 1);
                        fprintf(fid, '%u %f %f %f \n', dat');
                        fprintf(fid, 'End Coordinates\n');
                        
                        fprintf(fid, 'Elements\n');
                        ids = ce.getId();
                        nodeIds = arrayfun(@(x) x.getId(), nodes);
                        dat = sortrows([ids' nodeIds], 1);
                        format = ['%u ' repmat('%u ', 1, size(dat,2)-1) '\n'];
                        fprintf(fid, format, dat');
                        fprintf(fid, 'End Elements\n');
                    end
                    
                    %write model part nodes
                    nodes = mp.getNodes();
                    if ~isempty(nodes)
                        header = ['MESH "', mps{ii}, '_Nodes " dimension 3', ...
                            ' ElemType Point Nnode 1\n'];
                        fprintf(fid, header);
                        fprintf(fid, 'Coordinates\n');
                        ids = nodes.getId();
                        coords = nodes.getCoords();
                        dat = sortrows([ids' coords] ,1);
                        fprintf(fid, '%u %f %f %f \n', dat');
                        fprintf(fid, 'End Coordinates\n');
                        
                        %dummy nodal elements
                        fprintf(fid, 'Elements\n');
                        elids = nelms+1:nelms+length(ids);
                        dat = [elids; ids];
                        fprintf(fid, '%u %u\n', dat);
                        fprintf(fid, 'End Elements\n');
                        nelms = nelms+1+length(ids);
                    end
                end
            end
            
            fclose(fid);
        end
        
        function writeResults(obj, model)
            fid = fopen([obj.file, '.post.res'],'w');
            fprintf(fid, 'GiD Post Results File 1.0\n');
            nodeIds = model.getAllNodes().getId();
            freqs = model.getProcessInfo().getValue('FREQUENCY');
            if length(freqs) > 1
                steps = freqs;
            else
                steps = model.getProcessInfo().getValue('TIME');
            end
            
            elms = model.getAllElements();
            eTypes = string(arrayfun(@(x) class(x), elms, 'UniformOutput', 0));
            ueTypes = string(unique(eTypes));
            if ~isempty(obj.settings.elementalResults)
                for ii=1:length(ueTypes)
                    %identical element formulations with a different number
                    %of Gauss points are not supported here!
                    ce = elms(eTypes == ueTypes(ii));
                    cName = class(ce);
                    if isa(ce, 'TriangularElement')
                        [~, sep] = ce(1).getStressEvaluationPoints();
                    else
                        sep = ce(1).getStressEvaluationPoints();
                    end
                    
                    fprintf(fid, ['GaussPoints "', cName, '_gp" ElemType ', ...
                        ce.getElementType, '\n']);
                    fprintf(fid, ['Number of Gauss Points: ', ...
                        num2str(size(sep,1)), '\n']);
                    fprintf(fid, 'Natural Coordinates: Given\n');
                    fprintf(fid, '%f %f\n', sep');
                    fprintf(fid, 'End GaussPoints\n');
                end
            end
            
            for ii = 1:length(steps)
                for jj = 1:length(obj.settings.nodalResults)
                    name = char(obj.settings.nodalResults(jj));
                    type = checkPropertyName(name);
                    if strcmp(type, 'dofvec')
                        x = model.getAllNodes().getDofValue([name '_X'],ii);
                        y = model.getAllNodes().getDofValue([name '_Y'],ii);
                        try
                            z = model.getAllNodes().getDofValue([name '_Z'],ii);
                        catch
                            z = zeros(length(x),1);
                        end
                        if isreal(x) && isreal(y) && isreal(z)
                            resType = 'Vector';
                            dat = [nodeIds' x y z];
                            format ='%u %e %e %e\n';
                        else
                            resType = 'ComplexVector';
                            dat = [nodeIds' real(x) imag(x) real(y) imag(y) real(z) imag(z)];
                            format ='%u %e %e %e %e %e %e\n';
                        end
                    elseif strcmp(type, 'dofscalar')
                        val = model.getAllNodes().getDofValue(name,ii);
                        if isreal(val)
                            resType = 'Scalar';
                            dat = [nodeIds' val];
                            format = '%u %e\n';
                        else
                            resType = 'ComplexScalar';
                            dat = [nodeIds' real(val) imag(val)];
                            format = '%u %e %e\n';
                        end
                    elseif strcmp(type, 'variable1d')
                        val = model.getAllNodes().getPropertyValue(name);
                        if isreal(val)
                            resType = 'Scalar';
                            dat = [nodeIds' val];
                            format = '%u %e\n';
                        else
                            dat = [nodeIds' real(val) imag(val)];
                            format = '%u %e %e\n';
                            resType = 'ComplexScalar';
                        end
                    elseif strcmp(type, 'variable3d')
                        val = model.getAllNodes().getPropertyValue(name);
                        if isreal(val)
                            resType = 'Vector';
                            dat = [nodeIds' val];
                            format = '%u %e %e %e\n';
                        else
                            resType = 'ComplexScalar';
                            dat = [nodeIds' real(val(:,1)) imag(val(:,1)) ...
                                real(val(:,2)) imag(val(:,2)) ...
                                real(val(:,3)) imag(val(:,3))];
                        end
                    else
                        warning('unsupported result')
                    end
                    header = ['Result "' char(obj.settings.nodalResults(jj)) ...
                        '" "bm-mfem" ' num2str(steps(ii)) ' ', resType , ' OnNodes\n'];
                    fprintf(fid, header);
                    fprintf(fid, 'Values\n');
                    fprintf(fid, format, dat');
                    fprintf(fid, 'End Values\n');
                end %end nodal results
                
                
                for jj = 1:length(obj.settings.elementalResults)
                    if ii>1; warning('Stepped elemental results not supported yet\n'); end
                    name = char(obj.settings.elementalResults(jj));
                    type = checkPropertyName(name);
                    if strcmp(type, 'matrix')
                        for kk = 1:length(ueTypes)
                            ce = elms(eTypes == ueTypes(kk));
                            header = ['Result "' char(obj.settings.elementalResults(jj)) ...
                                '" "bm-mfem" ' num2str(steps(ii)) ...
                                ' Matrix OnGaussPoints "' class(ce) '_gp"\n'];
                            fprintf(fid, header);
                            fprintf(fid, 'Values\n');
                            for e = ce
                                fprintf(fid, num2str(e.getId));
                                %eDat = e.getPropertyValue(name,ii);
                                eDat = e.getPropertyValue(name);
                                dat = [squeeze(eDat(1,1,:)) squeeze(eDat(2,2,:)) ...
                                    squeeze(eDat(3,3,:)) squeeze(eDat(1,2,:)) ...
                                    squeeze(eDat(1,3,:)) squeeze(eDat(2,3,:))];
                                fprintf(fid, [repmat(' %e',1,6) '\n'], dat');
                            end
                            fprintf(fid, 'End Values\n');
                        end
                    else
                        warning('unsupported result')
                    end
                end
            end
            
            fclose(fid);
        end
        
    end
    
end

