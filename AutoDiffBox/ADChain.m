
classdef ADChain < handle
    % Chain rule of automatic differentiation
    
    properties
        target;
        factor;
        source;
        
        next;
    end
    
    
    methods(Static)
        
        function out = Top(in)
            persistent Top;
            if nargin > 0
                Top = in;
            end
            out = Top;
        end
        
        function [] = Clear()
            ADChain.Top([]);
        end
        
        function [] = Propagate(c)
            
            if (nargin == 0)
                c = ADChain.Top();
            end
            
            if ~isempty(c)
                c.source.adj = 1.0;
                
                while ~isempty(c)
                    % reverse accumulation
                    c.target.adj = c.target.adj + c.factor * c.source.adj;
                    
                    c = c.next;
                end
            end
            
            % ADChain.Clear();
        end
        
    end
    
    
    methods
        
        function obj = ADChain(f,s,t)
            obj.target = t;
            obj.factor = f;
            obj.source = s;
            
            obj.next = ADChain.Top();
            ADChain.Top(obj);
        end
        
        
    end
end


