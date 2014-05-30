function isIn = isMouseOver(element)
% ISMOUSEOVER checks if the mouse is over the element.
    mousePos = get( ...
        Gui.getParentFigure(element), ...
        'CurrentPoint' ...
    );

    elementPos = Gui.getGlobalPosition(element);
    
    isIn = ...
        mousePos(1) >= elementPos(1) && ...
        mousePos(1) <= (elementPos(1) + elementPos(3)) && ...
        mousePos(2) >= elementPos(2) && ...
        mousePos(2) <= (elementPos(2) + elementPos(4));
end