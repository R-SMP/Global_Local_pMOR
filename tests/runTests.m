function runTests( tests )
%RUNTESTS Summary of this function goes here
%   Detailed explanation goes here
if nargin ~= 1
    tests = 'small';
end

switch tests
    case 'small'
        run(AssemblerTests);
        run(ConditionTests);
        run(ElementTests);
        run(IOTests);
        run(SolverTests);
        run(VariousTests);
        
    case 'validation'
        run(ValidationTests);
        
    case 'all'
        run(ElementTests);
        run(ConditionTests);
        run(AssemblerTests);
        run(IOTests);
        run(SolverTests);
        run(VariousTests);
        run(ValidationTests);
        
    otherwise
        disp('possible input parameters are SMALL, VALIDATION, or ALL')

end

