x = 10.^(-3:0.01:1);
y = zeros(size(x));

parfor i = 1:numel(x)
    windowSize = max(300, ceil(x(i) * 300));
    T = windowSize / x(i);
    t = 0:(max(windowSize * 300, T * 300));
    y(i) = max( ...
        abs( ...
            Filter.median1d( ...
                sin(2*pi*t/T), ...
                windowSize ...
            ) ...
        ) ...
    );
end

ax = Paper.Axes();
ax.plot(x, y);
ax.setLogX
ax.setLogY();
ax.ylabel('Filter response');
ax.xlabel('Filter time \cdot frequency');
