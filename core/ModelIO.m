classdef ModelIO < handle
    %MODELIO Base class for all external input and output functionalities
    %   Detailed explanation goes here
    
    properties (Access = protected)
        file
        settings
    end
    
    methods
        function obj = ModelIO(file, settings)
            if nargin > 0
                obj.file = file;
                if nargin == 2
                    obj.settings = settings;
                else
                    obj.settings = struct();
                end
            end
        end
    end
    
    methods (Sealed)
        function setPreference(obj, setting, value)
            if isfield(obj.settings, setting)
                obj.settings.(setting) = value;
            else
                msg = [class(obj), ': The setting \"', ...
                    setting, '\" is not defined'];
                e = MException('MATLAB:bm_mfem:undefinedSetting',msg);
                throw(e);
            end
        end
        
        function checkIfInputExists(obj)
            if ~ exist(obj.file, 'file')
                msg = [class(obj), ': File ', obj.file, ' not found.'];
                e = MException('MATLAB:bm_mfem:fileNotFound',msg);
                throw(e);
            end
        end
    end
    
end
