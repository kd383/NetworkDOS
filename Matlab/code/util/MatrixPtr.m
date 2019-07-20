% Reference to a matrix variable.
classdef MatrixPtr < handle
    properties
        data
    end
    
    methods
        function obj = MatrixPtr(A)
            obj.data = A;
        end
        
        function d = diag(obj)
            d = diag(obj.data);
        end
        
        function t = trace(obj)
            t = trace(obj.data);
        end
    end
end