function newLSS_correlation(images, spmFile, rois, settings)
% FORMAT newLSS_correlation(images, spmFile, rois, fconnType, overwrite)
% Takes cell array of 4D images and, for each image, extracts mean beta
% series from mask and calculates voxel-wise correlation with that beta
% series. It then converts that R image to a Z image before taking the mean
% of Z images across all sessions for a given condition and saving that
% image.
%
%
% images:           Cell array of 4D images from
%                   generate_spm_singletrial_newLSS in format
%                   images{conds}{sessImages}
% spmFile:          Path to first level SPM.mat file from which LSS was
%                   performed. String.
% rois:             Cell array of paths to masks from which beta series
%                   will be extracted. For seed2voxel connectivity this
%                   beta series will be correlated with each of the voxels
%                   in the brain, while in roi2roi connectivity each roi's
%                   beta series will be correlated with every other roi's
%                   beta series. The latter is not set up yet.
% fconnType:        Which form of functional connectivity will be
%                   performed. Options are seed2voxel and roi2roi. String.

fconnType = settings.fconnType;
overwrite = settings.overwrite;

switch fconnType
    case 'seed2voxel'
        imageDir = fileparts(images{1}{1});
        parDir = fileparts(imageDir);
        outDir = [parDir '/correlations/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        
        load(spmFile);
        
        for iCond = 1:length(images)
            for jSess = 1:length(images{iCond})
                % Get rid of 4D from file name.
                [~, fileName] = fileparts(images{iCond}{jSess});
                fileName = fileName(3:end);
                for kROI = 1:length(rois)
                    % Get ROI name
                    [~, roiName] = fileparts(rois{kROI});
                    
                    % Create file names for condition outputs.
                    corrFilename{1} = [outDir '/Rcorr_' roiName fileName '.nii'];
                    corrFilename{2} = [outDir '/R_atanh_corr_' roiName fileName '.nii']; % Also known as z'
                    corrFilename{3} = [outDir '/Zcorr_' roiName fileName '.nii']; % Also known as z
                    zImages{iCond}{kROI}{jSess} = corrFilename{3};
                    
                    if overwrite || ~exist(corrFilename{3}, 'file')
                        % Correlate rest of brain with extracted timeseries from mask (R, R_atahn, Z).
                        correlation = beta_series_corr(images{iCond}{jSess}, SPM, rois{kROI});

                        % Write correlation (data) to corrFilename (name of output file).
                        write_corr_image(correlation{1}, corrFilename{1}, SPM.xVol);
                        write_corr_image(correlation{2}, corrFilename{2}, SPM.xVol);
                        write_corr_image(correlation{3}, corrFilename{3}, SPM.xVol);
                    else
                        fprintf('Exists: %s\n', corrFilename{3});
                    end
                end
            end
            for jROI = 1:length(rois)
                [outDir, fname] = fileparts(zImages{iCond}{jROI}{1});
                outName = [outDir '/mean_' fname(1:end-12) '.nii'];
                
                if overwrite || ~exist(outName, 'file')
                    create_summary_image(zImages{iCond}{jROI}, outName, 'mean(X)');
                else
                    fprintf('Exists: %s\n', outName);
                end
            end
        end
    case 'roi2roi'
        fprintf('This is not set up yet');
    otherwise
        error('I cannot make sense of this.');
end

end

%% Beta Series Correlation
function Cout = beta_series_corr(nii_loc, SPM, roi_loc, trimsd)
% FORMAT Cout = beta_series_corr(nii_loc, SPM, roi_loc, trimsd )
% Takes the beta series from one roi (the seed) and correlates it to
% all the voxels in the brain and saves the results as an image.
% Adapted from beta_series_correlation_nomars by Dennis Thompson
% without Events information.
% Calls spm_vol, roi_find_index, adjust_XYZ, and spm_get_data.
%
%
% nii_loc:          Path to LSS file (4D nifti) for particular condition. 
%                   Contains data to be correlated. String.
% SPM:              Loaded SPM variable from SPM.mat (first level from 
%                   original directory). Structure.
% roi_loc:          Path to mask files to be used. String.
% trimsd:           The number of standard deviations to use if you wish 
%                   the raw data to be Windsorized. Set to 0 if the raw 
%                   data are to be used. LSS default is 3. Double.

if ~exist('trimsd','var'), trimsd = 3; end
threshold = 0; % Find mask values greater than this

% Get header info for beta data.
V = spm_vol(nii_loc);

% Get ROI index and transform matrix.
if exist(roi_loc, 'file')
    [XYZ ROImat] = roi_find_index(roi_loc, threshold);
else
    error('Mask not found: %s\n', roi_loc);
end

% Generate XYZ locations for each beta image correcting for alignment
% issues and preallocate mean_roi vector.
betaXYZ = adjust_XYZ(XYZ, ROImat, V);
mean_roi = zeros(1, length(betaXYZ));

% Extract mean of ROI from each beta.
for iBeta = 1:length(betaXYZ),
    betasInROI = spm_get_data(V(iBeta), betaXYZ{iBeta});
    mean_roi(iBeta) = nanmean(betasInROI(:));
end

clear betasInROI

if trimsd > 0,
    mean_roi = trimts(mean_roi, trimsd);
end

% Get XYZ coordinates of all voxels from SPM.mat and then get beta
% values from LSS at those coordinates.
allDataAtBetas = spm_get_data(V, SPM.xVol.XYZ);
Cout = cell(3); Cout{1} = zeros(1, size(allDataAtBetas, 2));
for iVoxel = 1:size(allDataAtBetas, 2),
    if trimsd > 0,
        allDataAtBetas(:, iVoxel) = trimts(allDataAtBetas(:, iVoxel), trimsd);
    end
    Cout{1}(iVoxel) = corr(mean_roi', allDataAtBetas(:, iVoxel), 'type', 'Pearson');
end
Cout{2} = atanh(Cout{1});
stdZ = 1 / (length(V) - 3);
Cout{3} = Cout{2} / stdZ;
end

%% Write Correlation Image
function write_corr_image(Cout, imagename, xVol)
% FORMAT write_corr_image(Cout, imagename, xVol)
% Writes out whole brain correlation data into a nifti image. Verbatim
% Dennis Thompson's write_correlation_image.
%
%
% Cout:             A vector containing the correlations.
% imagename:        What the file is to be called, can include the full path.
% xVol:             A sub set of the SPM data structure.

nvoxels = size(xVol.XYZ, 2); % length of real data in brain - should match the length of Cout
corrData = NaN * ones(xVol.DIM(1), xVol.DIM(2), xVol.DIM(3)); % array of NaN size of brain volumne

% replace NaN's with corr results
for iVoxel = 1:nvoxels,
    corrData(xVol.XYZ(1, iVoxel), xVol.XYZ(2, iVoxel), xVol.XYZ(3, iVoxel)) = Cout(iVoxel);
end

[pathstr, name, ext] = fileparts(imagename);

% make nifti data structure see spm_vol - dt is [32 bit float, little endian]
CorrIm = struct ('fname', [name, ext], ...
    'dim',        [xVol.DIM'], ...
    'dt',         [16, 0], ...
    'mat',        xVol.M, ...
    'pinfo',      [1 0 0]', ...
    'descript',   'beta-correlation');

cwd = pwd;
cd(pathstr);

CorrIm = spm_create_vol(CorrIm, 'noopen');

[~] = spm_write_vol(CorrIm, corrData);

cd(cwd);
end

%% Trim stuff
function [y, ntrimmed] = trimts(y, sd)
% FORMAT [y, ntrimmed] = trimts(y, sd)
% Windsorizes data y by number of standard deviations sd.
% By Dennis Thompson.
%
% y:                1D vector of data.
% sd:               Number of standard deviations. Values out of this range 
%                   are to be replaced.
% y (output):       1D vectors of Windsorized data.
% ntrimmed:         Number of values replaced.

ntrimmed = 0;
idx = find(abs(y) > sd*std(y));
if ~isempty(idx),
    y(idx) = sign(y(idx)) * sd*std(y);
    ntrimmed = length(idx);
end
end

%% Extract Coordinates of ROI
function [index mat] = roi_find_index(ROI_loc, thresh)
% FORMAT [index mat] = roi_find_index(ROI_loc, thresh)
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
for n = 1:size(Y,3)
    % find values greater > thresh
    [xx yy] = find(squeeze(Y(:,:,n)) > thresh);
    if ~isempty(xx),
        zz = ones(size(xx))*n;
        index = [index,[xx';yy';zz']];
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

XYZ(4,:) = 1;
funcXYZ=cell(length(V));
for n = 1:length(V),
    if(iscell(V)),
       tmp = inv(V{n}.mat) * (ROImat * XYZ);
    else
        tmp = inv(V(n).mat) * (ROImat * XYZ);
    end
    funcXYZ{n} = tmp(1:3,:);
end
end