function h = plot(this, varargin)
    dataSizes = [this.dataSize];
    
    if (numel(this) == 1 || all(diff(dataSizes) == 0))
        time = [this.time];
        value = [this.value];
    else
        time = nan(max(dataSizes), numel(this));
        value = nan(max(dataSizes), numel(this));
        
        for i = 1:numel(this)
            time(1:dataSizes(i), i) = this(i).time;
            value(1:dataSizes(i), i) = this(i).value;
        end
    end
    
    if (any(size(time) == 0))
        time = nan(2, numel(this));
        value = nan(2, numel(this));
    end
    
    if (size(time, 1) == 1)
        time = vertcat(time, nan(1, numel(this)));
        value = vertcat(value, nan(1, numel(this)));
    end
    
    
    h = plot(time, value, 'DisplayName', {this.traceName}, varargin{:});
    
    for i = 1:numel(this)
%         ha = handle(h(i));
%         l = addlistener(this(i), 'change', @(~,~)update(ha, i));
%         addlistener(ha, 'ObjectBeingDestroyed', @(~,~)delete(l));
    end
    
    function update(h_, i)
        h_.YData = this(i).value;
        h_.XData = this(i).time;
        h_.DisplayName = this(i).traceName;
    end
end