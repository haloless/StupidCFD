
classdef ADScalar < handle
    % A single node in the AD chain
    
    properties
        val; % value
        adj; % adjoint
        
        chain;
    end
    
    methods
        
        %
        % constructor
        %
        function obj = ADScalar(v,a)
            if nargin <= 1
                a = 0;
            end
            obj.val = v;
            obj.adj = a;
            
            obj.chain = [];
        end
        
        %
        % operator overload
        %
        
        % -a
        function obj = uminus(obj1)
            v1 = double(obj1);
            v = -v1;
            obj = ADScalar(v);
            
            onUnaryOp(obj, obj1, -1.0);
        end
        
        % +a
        function obj = uplus(obj1)
            v1 = double(obj1);
            v = v1;
            obj = ADScalar(v);
            
            onUnaryOp(obj, obj1, +1.0);
        end
        
        % a + b
        function obj = plus(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = v1 + v2;
            obj = ADScalar(v);
            
            onBinaryOp(obj, obj1,1.0, obj2,1.0);
        end
        
        % a - b
        function obj = minus(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = v1 - v2;
            obj = ADScalar(v);
            
            onBinaryOp(obj, obj1,1.0, obj2,-1.0);
        end
        
        % a .* b
        function obj = times(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = v1 .* v2;
            obj = ADScalar(v);
            
            onBinaryOp(obj, obj1,v2, obj2,v1);
        end
        
        % a ./ b
        function obj = rdivide(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = v1 ./ v2;
            obj = ADScalar(v);
            
            coef1 = 1 ./ v2;
            coef2 = -v1 ./ v2.^2;
            onBinaryOp(obj, obj1,coef1, obj2,coef2);
        end
        
        % a .^ b
        function obj = power(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = v1 .^ v2;
            obj = ADScalar(v);
            
            coef1 = v2 .* v1.^(v2-1);
            coef2 = v .* log(v1);
            onBinaryOp(obj, obj1,coef1, obj2,coef2);
        end
        
        % a * b
        function obj = mtimes(obj1,obj2)
            obj = times(obj1,obj2);
        end
        
        % a / b
        function obj = mrdivide(obj1,obj2)
            obj = rdivide(obj1,obj2);
        end
        
        % a ^ b
        function obj = mpower(obj1,obj2)
            obj = power(obj1,obj2);
        end
        
        %
        function ret = lt(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            ret = (v1 < v2);
        end
        
        function ret = gt(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            ret = (v1 > v2);
        end
        
        %
        % function overload
        %
        function obj = sin(obj1)
            v1 = double(obj1);
            v = sin(v1);
            obj = ADScalar(v);
            
            a = cos(v1);
            onUnaryOp(obj, obj1, a);
        end
        
        function obj = cos(obj1)
            v1 = double(obj1);
            v = cos(v1);
            obj = ADScalar(v);
            
            a = -sin(v1);
            onUnaryOp(obj, obj1, a);
        end
        
        function obj = tan(obj1)
            v1 = double(obj1);
            v = tan(v1);
            obj = ADScalar(v);
            
            a = 1 ./ cos(v1).^2;
            onUnaryOp(obj, obj1, a);
        end
        
        % atan(x)
        function obj = atan(obj1)
            v1 = double(obj1);
            v = atan(v1);
            obj = ADScalar(v);
            
            a = 1 ./ (v.^2 + 1);
            onUnaryOp(obj, obj1, a);
        end
        
        % atan2(y,x)
        function obj = atan2(obj1,obj2)
            v1 = double(obj1);
            v2 = double(obj2);
            v = atan2(v1,v2);
            obj = ADScalar(v);
            
            a1 = +v2 ./ (v1.^2 + v2.^2);
            a2 = -v1 ./ (v1.^2 + v2.^2);
            onBinaryOp(obj, obj1,a1, obj2,a2);
        end
        
        %
        % helper routines
        %
        function d = double(obj)
            d = obj.val;
        end
        
        function [] = onBinaryOp(obj, lvar,lcoef, rvar,rcoef)
            if isa(lvar,'ADScalar')
                lop = ADChain(lcoef, obj, lvar);
            end
            if isa(rvar,'ADScalar')
                rop = ADChain(rcoef, obj, rvar);
            end
            
            obj.chain = ADChain.Top();
        end
        
        function [] = onUnaryOp(obj, var, coef)
            if isa(var,'ADScalar')
                op = ADChain(coef, obj, var);
            end
            
            obj.chain = ADChain.Top();
        end
        
        %
        function [] = propagate(obj)
            c = obj.chain;
            ADChain.Propagate(c);
        end
        
    end
    
end % ADScalar









