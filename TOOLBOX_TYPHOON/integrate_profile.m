function [ i ] = integrate_profile( x, y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    

    close all

    %xlim = [x(1) x(end)]; 
    xlim = [1 length(y)]; 
    plot( xlim(1):xlim(2), y, 'b'), hold on

    ylim = [min(y) max(y)];
    set(gca, 'XLim', xlim)
    set(gca, 'YLim', ylim)
    title('Double-click on the red line to finish')
    h1 = imline(gca, [1 1]*((xlim(2)-xlim(1))*0.2+xlim(1)), ylim);
    h2 = imline(gca, [1 1]*((xlim(2)-xlim(1))*0.8+xlim(1)), ylim);

    setColor(h1,[1 0 0]);
    setColor(h2,[0 1 0]);

    setPositionConstraintFcn(h1, @(pos)[ [pos(2,1); pos(2,1)] ylim'  ])
    setPositionConstraintFcn(h2, @(pos)[ [pos(2,1); pos(2,1)] ylim'  ])
    
    pos_line1 = wait(h1);
    pos_line2 = getPosition(h2);

    i1 = round(pos_line1(1,1));
    i2 = round(pos_line2(1,1));
    
    i = [min([i1 i2]) max([i1 i2])];
    close all


end

