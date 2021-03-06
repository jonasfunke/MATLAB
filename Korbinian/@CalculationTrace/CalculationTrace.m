classdef CalculationTrace < AbstractTrace & handle
    %CalculationTrace
    
    properties(SetAccess=private)
        trace1
        operation
        trace2
    end
    
    properties(Access=private,Transient)
        calculatedValue
    end
    
    methods
        function this = CalculationTrace(trace1, operation, trace2)
            
            this.trace1 = trace1;
            this.operation = operation;
            if ( ...
                numel(trace1.time) ~= numel(trace2.time) || ...
                any(trace1.time ~= trace2.time) ...
            )
                trace2 = ResampledTrace(trace2);
                trace2.setResampledTime(trace1.time);
            end
            this.trace2 = trace2;
            
            this.registerListeners();
        end
        
        function registerListeners(this)
            
            
            l = [ ...
                addlistener(this.trace1, 'change', @this.resetCalculated), ...
                addlistener(this.trace2, 'change', @this.resetCalculated) ...
            ];
            
            addlistener(this, 'ObjectBeingDestroyed', @(~,~)delete(l));
        end
        
        function time = getTime(this)
            time = this.trace1.time;
        end
        
        function value = getValue(this)
            if (isempty(this.calculatedValue))
                this.calculate();
            end
            value = this.calculatedValue;
        end
        
        function timeUnit = getTimeUnit(this)
            timeUnit = this.trace1.getTimeUnit();
        end
        function timeName = getTimeName(this)
            timeName = this.trace1.getTimeName();
        end
        function valueUnit = getValueUnit(this)
            valueUnit = this.trace1.getValueUnit();
        end
        function valueName = getValueName(this)
            valueName = this.trace1.getValueName();
        end
        function name = getName(this)
            name = this.trace1.getName();
        end
    end
    
    methods (Access=private)
        function resetCalculated(this, ~, ~)
            for o = this
                o.calculatedValue = [];
                o.notify('change');
            end
        end
        function calculate(this)
            for o = this
                value1 = o.trace1.value;
                value2 = o.trace2.value;
                switch o.operation
                    case '+'
                        value = value1 + value2;
                    case '-'
                        value = value1 - value2;
                    case '*'
                        value = value1 .* value2;
                    case '/'
                        value = value1 ./ value2;
                    case '^'
                        value = value1 .^ value2;
                end
                
                o.calculatedValue = value;
            end
        end
    end
end

