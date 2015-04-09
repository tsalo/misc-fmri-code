function [outString, baString] = getBaAndRegion(input)
% FORMAT [outString, baString] = getBaAndRegion(clusterImage)
% Takes cluster matrix, determines regions.
spmDir = fileparts(which('spm'));
spmDirs = strsplit(spmDir, '/');
spmDirs = spmDirs(1:end-1);
spmDir = strjoin(spmDirs, '/');
spmDir = fullfile(spmDir, 'spm8');
pickatlasDir = fullfile(spmDir, 'toolbox/wfu_pickatlas');

typeNiiFile = fullfile(pickatlasDir, 'MNI_atlas_templates/TD_type.nii');
typeTextFile = fullfile(pickatlasDir, 'MNI_atlas_templates/TD_type.txt');
baNiiFile = fullfile(pickatlasDir, 'MNI_atlas_templates/TD_brodmann.nii');
baTextFile = fullfile(pickatlasDir, 'MNI_atlas_templates/TD_brodmann.txt');
labelNiiFile = fullfile(pickatlasDir, 'MNI_atlas_templates/atlas116.nii');
labelTextFile = fullfile(pickatlasDir, 'MNI_atlas_templates/atlas116.txt');

[typeHeader, typeValues] = loadNifti(typeNiiFile);
tissueTypes = loadTableFile(typeTextFile);

[~, baValues] = loadNifti(baNiiFile);
baTypes = loadTableFile(baTextFile);

[~, labelValues] = loadNifti(labelNiiFile);
labelTypes = loadTableFile(labelTextFile);

if isa(input, 'string')
    clusterImage = input;
    [xyz, maskMatrix] = roiFindIndex(clusterImage, 0);
    clusterXyz = adjustXyz(xyz, maskMatrix, typeHeader);
    clusterIndex = sub2ind(size(typeValues), clusterXyz{1}(1, :), clusterXyz{1}(2, :), clusterXyz{1}(3, :))';

    clusterTypes = typeValues(clusterIndex);
    clusterBas = baValues(clusterIndex);
    clusterLabels = labelValues(clusterIndex);
elseif isa(input, 'double')
    coordinates = input;
    [xMm, yMm, zMm] = vox2mm(coordinates, typeNiiFile);
    
    clusterTypes = typeValues(xMm, yMm, zMm);
    clusterBas = baValues(xMm, yMm, zMm);
    clusterLabels = labelValues(xMm, yMm, zMm);
end

% Tissue Types
result = setdiff(unique(clusterTypes), 0);
instances = histc(clusterTypes, result);
[~, index] = sort(instances, 'descend');
sortedValues = result(index);
typeCells = cell(size(sortedValues))';
for iType = 1:length(sortedValues)
    ttIndex = strcmp(tissueTypes(:, 1), num2str(sortedValues(iType)));
    tissueType = tissueTypes{ttIndex, 2};
    typeCells{iType} = tissueType;
end
if ~isempty(typeCells)
    typeString = strjoin(typeCells, '| ');
else
    typeString = 'Undefined';
end

% Labels
values = setdiff(unique(clusterLabels), 0);
instances = histc(clusterLabels, values);
[~, index] = sort(instances, 'descend');
sortedValues = values(index);
labelCells = cell(size(sortedValues))';
for iLabel = 1:length(sortedValues)
    labelIndex = strcmp(labelTypes(:, 1), num2str(sortedValues(iLabel)));
    labelType = labelTypes{labelIndex, 2};
    labelCells{iLabel} = labelType;
end
if ~isempty(labelCells)
    labelString = strjoin(labelCells, '| ');
else
    labelString = 'Undefined';
end

% BAs
values = setdiff(unique(clusterBas), 0);
instances = histc(clusterBas, values);
[~, index] = sort(instances, 'descend');
sortedValues = values(index);
baCells = cell(size(sortedValues))';
for iBa = 1:length(sortedValues)
    baIndex = strcmp(baTypes(:, 1), num2str(sortedValues(iBa)));
    baType = strrep(baTypes{baIndex, 2}, 'brodmann area ', '');
    baCells{iBa} = baType;
end
if ~isempty(baCells)
    baString = strjoin(baCells, '| ');
else
    baString = 'Undefined';
end

if ~strfind(typeString, 'Gray')
    outString = typeString;
else
    outString = labelString;
end
end

function [xSub, ySub, zSub] = vox2mm(varargin)
% FORMAT [xSub, ySub, zSub] = vox2mm(x, y, z, niiFile)
% or
% FORMAT [xSub, ySub, zSub] = vox2mm([x, y, z], niiFile)
% Converts voxel coordinates to subscript values denoting position in data
% matrix.
%
% Inputs:
% x:        X-coordinate.
% y:        Y-coordinate.
% z:        Z-coordinate.
% niiFile:  Nifti file from which alignment matrix will be used.
%
% Outputs:
% xSub:     X subscript.
% ySub:     Y subscript.
% zSub:     Z subscript.
%
% 150305 Created by Taylor Salo

if nargin == 2
    loc = varargin{1};
    niiFile = varargin{2};
elseif nargin == 4
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    loc = [x, y, z];
    niiFile = varargin{4};
else
    error('You gave bad inputs.');
end

v2m = spm_get_space(niiFile);
mm = (loc - v2m(1:3, 4)') / v2m(1:3, 1:3);
xSub = mm(1);
ySub = mm(2);
zSub = mm(3);
end
