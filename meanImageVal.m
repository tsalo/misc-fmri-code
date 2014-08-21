function meanVal = meanImageVal(valuesFile, maskFile)
% FORMAT meanVal = meanImageVal(valuesFile, maskFile)
% Gives you the mean across a ResMS file.

Vvals = spm_vol(valuesFile);

if exist('maskFile', 'var')
    Vmask = spm_vol(maskFile);
    [maskRawXYZ maskMat] = roi_find_index(maskFile, 0);
    maskXYZ = adjust_XYZ(maskRawXYZ, maskMat, Vvals);
    imageVals = spm_get_data(valuesFile, maskXYZ{1});
    imageVals(isnan(imageVals)) = 0;
    meanVal = mean(imageVals);
else
    disp('Warning. No mask given. Extracting mean from whole volume.')
    imageVals = spm_read_vols(Vvals);
    imageVals(isnan(imageVals)) = 0;
    meanVal = mean(imageVals(:));
end
end

%% Extract Coordinates of ROI
function [index, mat] = roi_find_index(ROI_loc, thresh)
% FORMAT [index, mat] = roi_find_index(ROI_loc, thresh)
% Returns the XYZ address of voxels with values greater than threshold. 
% By Dennis Thompson.
% 
%
% ROI_loc:          String pointing to nifti image.
% thresh:           Threshold value, defaults to zero. Double.
if ~exist('thresh','var'),
    thresh = 0;
end

data = nifti(ROI_loc);
Y = double(data.dat);
Y(isnan(Y)) = 0;
index = [];
for n = 1:size(Y, 3)
    % find values greater > thresh
    [xx, yy] = find(squeeze(Y(:, :, n)) > thresh);
    if ~isempty(xx),
        zz = ones(size(xx)) * n;
        index = [index, [xx'; yy'; zz']];
    end
end

mat = data.mat;
end

%% Adjust Coordinates of ROI
function funcXYZ = adjust_XYZ(XYZ, ROImat, V)
% FORMAT funcXYZ = adjust_XYZ(XYZ, ROImat, V)
% By Dennis Thompson.
%
%
% XYZ:              Output from lss_roi_find_index.
% ROImat:           Output from lss_roi_find_index.
% V:                Header information of nifti file from spm_vol.
XYZ(4, :) = 1;
funcXYZ = cell(length(V));
for n = 1:length(V),
    if(iscell(V)),
       tmp = inv(V{n}.mat) * (ROImat * XYZ);
    else
        tmp = inv(V(n).mat) * (ROImat * XYZ);
    end
    funcXYZ{n} = tmp(1:3, :);
end
end

