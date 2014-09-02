function newLSS_correlation(images, spmFile, rois, settings)
% FORMAT newLSS_correlation(images, spmFile, rois, settings)
% Takes cell array of 4D images and, for each image, extracts mean beta
% series from masks and computes one of two forms of functional
% connectivity: seed2voxel or roi2roi. In seed2voxel functional
% connectivity, the function calculates voxel-wise correlation between
% mean beta series from a given mask and  beta series from the rest of the
% brain. It then converts that R image to a Z image before taking the mean
% of Z images across all sessions for a given condition and saving that
% image. In roi2roi functional connectivity, it extracts the mean beta
% series from each mask given and correlates those beta series to produce a
% matrix, which it then averages across sessions. The outputted matrix is
% also in Z values. I have not done much research into graph theory so
% nothing is included for analyzing the roi2roi results.
%
%
% images:             Cell array of 4D images from
%                     generate_spm_singletrial_newLSS in format
%                     images{conds}{sessImages}
% spmFile:            Path to first level SPM.mat file from which LSS was
%                     performed. String.
% rois:               Cell array of paths to masks from which beta series
%                     will be extracted. For seed2voxel connectivity this
%                     beta series will be correlated with each of the 
%                     voxels in the brain, while in roi2roi connectivity
%                     each roi's beta series will be correlated with every
%                     other roi's beta series.
% settings:           Additional settings. Structure.
% settings.fConnType: Which form of functional connectivity will be
%                     performed. Options are "seed2voxel" and "roi2roi".
%                     String.
% settings.overwrite: Overwrite any pre-existing files (1) or not (0).
%                     Double.
switch settings.fConnType
    case 'seed2voxel'
        imageDir = fileparts(images{1}{1});
        parDir = fileparts(imageDir);
        outDir = [parDir '/seed2voxel_correlations/'];
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
                    
                    if settings.overwrite || ~exist(corrFilename{3}, 'file')
                        % Correlate rest of brain with extracted timeseries from mask (R, R_atahn, Z).
                        meanROI = extract_beta_series(images{iCond}{jSess}, rois{kROI});
                        correlation = beta_series_corr(images{iCond}{jSess}, SPM, meanROI);

                        % Write correlation (data) to corrFilename (name of output file).
                        write_corr_image(correlation{1}, corrFilename{1}, SPM.xVol);
                        write_corr_image(correlation{3}, corrFilename{3}, SPM.xVol);
                    else
                        fprintf('Exists: %s\n', corrFilename{3});
                    end
                end
            end
            for jROI = 1:length(rois)
                [outDir, fname] = fileparts(zImages{iCond}{jROI}{1});
                outName = [outDir '/mean_' fname(1:end-8) '.nii'];
                if settings.overwrite || ~exist(outName, 'file')
                    % If there are multiple runs (i.e. LSS), create a mean file.
                    if length(zImages{iCond}{jROI}) > 1
                        create_summary_image(zImages{iCond}{jROI}, outName, 'mean(X)');
                    end
                else
                    fprintf('Exists: %s\n', outName);
                end
            end
        end
    case 'roi2roi'
        imageDir = fileparts(images{1}{1});
        parDir = fileparts(imageDir);
        outDir = [parDir '/roi2roi_correlations/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        
        load(spmFile);
        
        for iCond = 1:length(images)
            for jSess = 1:length(images{iCond})
                % Get rid of 4D from file name.
                [~, fileName, ~] = fileparts(images{iCond}{jSess});
                fileName = fileName(3:end);
                
                % Create file names for condition outputs.
                corrMatrixName = [outDir '/Zcorr' fileName '.mat'];
                zImages{iCond}{jSess} = corrMatrixName;

                if settings.overwrite || ~exist(corrMatrixName, 'file')
                    zCorrMatrix = zeros(length(rois));
                    roiNames = cell(length(rois), 1);
                    roiBetaSeries = cell(length(rois), 1);
                    for kROI = 1:length(rois)
                        % Get ROI name
                        [~, roiNames{kROI}] = fileparts(rois{kROI});
                        roiBetaSeries{kROI} = extract_beta_series(images{iCond}{jSess}, rois{kROI});
                    end
                    nTrials = length(roiBetaSeries{1});
                    for kROI = 1:length(rois)
                        for mROI = 1:length(rois)
                            rCorr = corr(roiBetaSeries{kROI}', roiBetaSeries{mROI}', 'type', 'Pearson');
                            stdZ = 1 / (nTrials - 3);
                            zCorrMatrix(kROI, mROI) = atanh(rCorr) / stdZ;
                        end
                    end
                    corrStruct.rois = roiNames;
                    corrStruct.corrMatrix = zCorrMatrix;
                    save(corrMatrixName, 'corrStruct'); 
                end
            end
            [outDir, fname] = fileparts(zImages{iCond}{1});
            outName = [outDir '/mean_' fname(1:end-8) '.mat'];
            
            if settings.overwrite || ~exist(outName, 'file')
                % If there are multiple runs (i.e. LSS), create a mean file.
                if length(zImages{iCond}) > 1
                    allCorr = [];
                    for jSess = 1:length(images{iCond})
                        load(zImages{iCond}{jSess});
                        allCorr(:, :, jSess) = corrStruct.corrMatrix;
                        meanStruct.rois = corrStruct.rois;
                        clear corrStruct
                    end
                    meanStruct.meanCorr = mean(allCorr, 3);
                    save(outName, 'meanStruct');
                end
            else
                fprintf('Exists: %s\n', outName);
            end
        end
    otherwise
        error(['I cannot make sense of this: ' settings.fConnType]);
end
end

%% Beta Series Correlation
function meanROI = extract_beta_series(niiLoc, roiLoc, trimsd)
% FORMAT meanROI = extract_beta_series(niiLoc, roiLoc, trimsd)
% Extracts and returns mean beta series from 4D niiLoc within ROI.
% Adapted from beta_series_correlation_nomars by Dennis Thompson
% without Events information.
% Calls spm_vol, roi_find_index, adjust_XYZ, and spm_get_data.
%
%
% niiLoc:           Path to LSS file (4D nifti) for particular condition. 
%                   Contains data to be correlated. String.
% roiLoc:           Path to mask file to be used. String.
% trimsd:           The number of standard deviations to use if you wish 
%                   the raw data to be Windsorized. Set to 0 if the raw 
%                   data are to be used. LSS default is 3. Double.

if ~exist('trimsd', 'var'), trimsd = 3; end
threshold = 0; % Find mask values greater than this

% Get header info for beta data.
V = spm_vol(niiLoc);

% Get ROI index and transform matrix.
if exist(roiLoc, 'file')
    [XYZ ROImat] = roi_find_index(roiLoc, threshold);
else
    error('Mask not found: %s\n', roiLoc);
end

% Generate XYZ locations for each beta image correcting for alignment
% issues and preallocate meanROI vector.
betaXYZ = adjust_XYZ(XYZ, ROImat, V);
meanROI = zeros(1, length(betaXYZ));

% Extract mean of ROI from each beta.
for iBeta = 1:length(betaXYZ),
    betasInROI = spm_get_data(V(iBeta), betaXYZ{iBeta});
    meanROI(iBeta) = nanmean(betasInROI(:));
end

clear betasInROI

if trimsd > 0,
    meanROI = trimts(meanROI, trimsd);
end
end

%% Create Correlation image
function Cout = beta_series_corr(niiLoc, SPM, meanROI, trimsd)
% FORMAT Cout = beta_series_corr(niiLoc, SPM, meanROI, trimsd)
% Correlates meanROI beta series with beta series of all other voxels in
% brain. Converts resulting R matrix to Z matrix.
% Adapted from beta_series_correlation_nomars by Dennis Thompson
% without Events information.
% Calls spm_vol and spm_get_data.
%
%
% niiLoc:           Path to LSS file (4D nifti) for particular condition. 
%                   Contains data to be correlated. String.
% SPM:              Loaded SPM variable from SPM.mat (first level from 
%                   original directory). Structure.
% meanROI:          Mean beta series correlated with all other voxels' beta
%                   series in brain. Double vector.
% trimsd:           The number of standard deviations to use if you wish 
%                   the raw data to be Windsorized. Set to 0 if the raw 
%                   data are to be used. LSS default is 3. Double.
if ~exist('trimsd', 'var'), trimsd = 3; end
V = spm_vol(niiLoc);
allDataAtBetas = spm_get_data(V, SPM.xVol.XYZ);
Cout = cell(3); Cout{1} = zeros(1, size(allDataAtBetas, 2));
for iVoxel = 1:size(allDataAtBetas, 2),
    if trimsd > 0,
        allDataAtBetas(:, iVoxel) = trimts(allDataAtBetas(:, iVoxel), trimsd);
    end
    Cout{1}(iVoxel) = corr(meanROI', allDataAtBetas(:, iVoxel), 'type', 'Pearson');
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

%% Trim Timeseries
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
funcXYZ = cell(length(V));
for n = 1:length(V),
    if(iscell(V)),
       tmp = inv(V{n}.mat) * (ROImat * XYZ);
    else
        tmp = inv(V(n).mat) * (ROImat * XYZ);
    end
    funcXYZ{n} = tmp(1:3,:);
end
end

%% Write Out Mean Image
function create_summary_image(cellVols, outName, expr)
% FORMAT create_summary_image(cellVols, outName, expr)
% Calls spm_imcalc to perform calculations (sum, mean, std) on a cell array
% of volumes to output one summary volume.
%
% cellVols:         Cell array of images upon which calculations will be
%                   performed.
% outName:          The name (with full path) of the image to be written. String.
% expr:             The calculation to be performed upon the images in cellVols.
%                   String. Possible values: sum(X), mean(X), std(X), maybe others.

for iVol = 1:length(cellVols)
    Vi(iVol) = spm_vol(cellVols{iVol});
end
Vo = Vi(1);
Vo.fname = outName;

spm_imcalc(Vi, Vo, expr, {1, 0, 0});
end
