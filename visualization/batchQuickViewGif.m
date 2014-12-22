% Make images!
inputFile = '/home/tsalo/fmri-images-rise/tcan_final.nii';
outFolder = '/home/tsalo/fmri-images-rise/';
view = 'sag';

switch view
    case 'ax'
        % zMin = -50
        % zMax = 84
        zCoords = [-50:2:84 82:-2:-48]';
        yCoords = zeros(size(zCoords));
        xCoords = zeros(size(zCoords));
        coordsVary = zCoords;
    case 'cor'
        % yMin = -112
        % yMax = 76
        yCoords = [-112:2:76 74:-2:-110]';
        xCoords = zeros(size(yCoords));
        zCoords = zeros(size(yCoords));
        coordsVary = yCoords;
    case 'sag'
        % xMin = -78
        % xMax = 78
        xCoords = [-78:2:78 76:-2:-76]';
        yCoords = zeros(size(xCoords));
        zCoords = zeros(size(xCoords));
        coordsVary = xCoords;
end

coordinates = [xCoords yCoords zCoords];

for iSlice = 1:length(coordsVary)
    saveFile = [outFolder sprintf('%03d', iSlice) '_' num2str(coordsVary(iSlice))];
    quickViewGif(inputFile, saveFile, coordinates(iSlice, :), view)
end