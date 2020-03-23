function output_txt = cursorOutputFunction(~, event_obj, elecLabels, labels, values)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

if strcmp(event_obj.Target.Type, 'scatter')
    idxElec2 = get(event_obj, 'DataIndex');
    if ~isempty(labels)
        output_txt = {elecLabels{idxElec2}, labels{idxElec2}};
    elseif ~isempty(values)
        output_txt = {elecLabels{idxElec2}, sprintf('%.3f', values(idxElec2))};
    else
        output_txt = {elecLabels{idxElec2}};
    end
else
    pos = get(event_obj, 'Position');
    output_txt = {['X: ',num2str(pos(1),4)], ['Y: ',num2str(pos(2),4)], ['Z: ',num2str(pos(3),4)]};
end
