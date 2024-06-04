function []=export_to_ppt()

% Create a PowerPoint presentation
ppt = exportToPPTX('new');

% Get handles to open figures
figHandles = findobj('Type', 'figure');

% Loop through each figure and export to PowerPoint
for i = 1:numel(figHandles)
    % Copy figure to clipboard
    figHandle = figHandles(i);
    figCopy = copyobj(figHandle, 0);
    
    % Add figure to PowerPoint slide
    ppt.addSlide('Layout', 'Title and Content');
    ppt.addText(['Figure ', num2str(i)], 'Position', 'title');
    ppt.addPicture(gcf, 'Position', 'content');
    
    % Close the copied figure
    delete(figCopy);
end

% Save the PowerPoint presentation
ppt.saveAs('MyFigures.pptx');

% Close the PowerPoint presentation
ppt.close();

end