function str = char(this)
    str = sprintf( ...
        'sliced: %.2f - %.2f', ...
        this.startTime, ...
        this.endTime ...
    );
end