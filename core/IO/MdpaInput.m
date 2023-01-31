classdef MdpaInput < ModelIO
    %MDPAINPUT Read model input from mdpa
    %   Detailed explanation goes here

    properties (Access = private)
    end

    methods

        function obj = MdpaInput(file)
            settings = struct('readConditions', true);

            if nargin == 0
                super_args = {};
            elseif nargin == 1
                super_args = {file, settings};
            end

            obj@ModelIO(super_args{:});
            obj.checkIfInputExists();
        end

        function model = readModel(obj)
            fid = fopen(obj.file);
            props = PropertyContainer;
            model = FemModel;

            tline = fgetl(fid);

            while ischar(tline)
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end

                if startsWith(tline,'Begin Properties')
                    nProp = strsplit(tline);
                    nProp = str2double(nProp{end}) + 1; %MDPA starts at 0
                    [fid,props] = obj.readProperties(fid,props,nProp);
                end

                if startsWith(tline,'Begin Nodes')
                    [fid,model] = obj.readNodes(fid,model);
                end

                if startsWith(tline,'Begin Elements')
                    etype = strsplit(tline);
                    etype = cell2mat(etype(3));
                    etype = erase(etype,'//');
                    [fid,model] = obj.readElements(fid,model,etype,props);
                end

                if startsWith(tline,'Begin Conditions') && obj.settings.readConditions
                    ctype = strsplit(tline);
                    ctype = cell2mat(ctype(3));
                    ctype = erase(ctype,'//');
                    [fid,model] = obj.readConditions(fid,model,ctype,props);
                end

                if startsWith(tline,'Begin SubModelPart')
                    name = strsplit(tline);
                    name = cell2mat(name(3));
                    [fid,model] = obj.readSubModelParts(fid,model,name,obj.settings.readConditions);
                end

                tline = fgetl(fid);
            end

        end

    end

    methods (Static)

        function [fid,props] = readProperties(fid,props,nProp)
            tline = fgetl(fid);
            property = PropertyContainer;
            while ~ startsWith(tline,'End Properties')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end

                propData = strsplit(tline);
                propData(cellfun('isempty',propData)) = [];
                property.addValue(cell2mat(propData(1)),str2double(propData(2)));
                tline = fgetl(fid);
            end
            props(nProp) = property;
        end

        function [fid,model] = readNodes(fid,model)
            dat = textscan(fid,'%u%f%f%f','delimiter','\t','MultipleDelimsAsOne',true);
            for ii = 1:length(dat{1})
                model.addNewNode(dat{1}(ii),dat{2}(ii),dat{3}(ii),dat{4}(ii));
            end
        end

        function [fid,model] = readElements(fid,model,etype,props)
            tline = fgetl(fid);
            while ~ startsWith(tline,'End Elements')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end

                e = cell2mat(textscan(tline,'%f'));
                model.addNewElement(etype,e(1),e(3:end),props(e(2)+1));
                tline = fgetl(fid);
            end
        end

        function [fid,model] = readConditions(fid,model,ctype,props)
            tline = fgetl(fid);
            while ~startsWith(tline,'End Conditions')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end
                c = cell2mat(textscan(tline,'%f'));
                model.addNewCondition(ctype,c(1),c(3:end),props(c(2)+1));
                tline = fgetl(fid);
            end
        end

        function [fid,model] = readSubModelParts(fid,model,name,readConditions)
            tline = fgetl(fid);
            nodes = [];
            elements = [];
            conditions = [];
            while ~ startsWith(tline,'End SubModelPart')
                if startsWith(strtrim(tline),'//'); tline = fgetl(fid); continue; end

                if contains(tline,'Begin SubModelPartNodes')
                    tline = fgetl(fid);
                    while ~ contains(tline,'End SubModelPartNodes')
                        nodes = [nodes str2double(tline)]; %#ok<AGROW>
                        tline = fgetl(fid);
                    end
                end

                if contains(tline,'Begin SubModelPartElements')
                    tline = fgetl(fid);
                    while ~ contains(tline,'End SubModelPartElements')
                        elements = [elements str2double(tline)]; %#ok<AGROW>
                        tline = fgetl(fid);
                    end
                end

                if contains(tline,'Begin SubModelPartConditions') && readConditions
                    tline = fgetl(fid);
                    while ~ contains(tline,'End SubModelPartConditions')
                        conditions = [conditions str2double(tline)]; %#ok<AGROW>
                        tline = fgetl(fid);
                    end
                end

                tline = fgetl(fid);
            end

            model.addNewModelPart(name,nodes,elements,conditions);
        end

    end

end

