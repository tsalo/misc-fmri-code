function newLSS_correlation(images, rois, settings)
% FORMAT newLSS_correlation(images, rois, settings)
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
                        meanRoi = extractBetaSeries(images{iCond}{jSess}, rois{kROI})';
                        
                        niiHeader = spm_vol(images{iCond}{jSess});
                        [Y, ~] = spm_read_vols(niiHeader(1));
                        nanIdx = find(~isnan(Y));
                        [x, y, z] = ind2sub(size(Y), nanIdx);
                        xVol.XYZ = [x y z].';
                        xVol.DIM = size(Y);
                        xVol.M = niiHeader.mat;
                        
                        correlation = correlateBetaSeries(niiHeader, xVol, meanRoi);

                        % Write correlation (data) to corrFilename (name of output file).
                        writeCorrelationImage(correlation{1}, corrFilename{1}, xVol);
                        writeCorrelationImage(correlation{3}, corrFilename{3}, xVol);
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
                        writeSummaryImage(zImages{iCond}{jROI}, outName, 'mean(X)');
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
        
        for iCond = 1:length(images)
            for jSess = 1:length(images{iCond})
                % Get rid of 4D from file name.
                [~, fileName, ~] = fileparts(images{iCond}{jSess});
                fileName = fileName(3:end);
                
                % Create file names for condition outputs.
                corrMatrixName{2} = [outDir '/Rcorr' fileName '.mat'];
                corrMatrixName{2} = [outDir '/Zcorr' fileName '.mat'];
                zImages{iCond}{jSess} = corrMatrixName;

                if settings.overwrite || ~exist(corrMatrixName, 'file')
                    rCorrMatrix = zeros(length(rois));
                    zCorrMatrix = zeros(length(rois));
                    roiNames = cell(length(rois), 1);
                    roiBetaSeries = cell(length(rois), 1);
                    for kROI = 1:length(rois)
                        % Get ROI name
                        [~, roiNames{kROI}] = fileparts(rois{kROI});
                        roiBetaSeries{kROI} = extractBetaSeries(images{iCond}{jSess}, rois{kROI})';
                    end
                    nTrials = length(roiBetaSeries{1});
                    for kROI = 1:length(rois)
                        for mROI = 1:length(rois)
                            rCorrMatrix(kROI, mROI) = corr(roiBetaSeries{kROI}, roiBetaSeries{mROI}, 'type', 'Pearson');
                            stdZ = 1 / (nTrials - 3);
                            zCorrMatrix(kROI, mROI) = atanh(rCorrMatrix(kROI, mROI)) / stdZ;
                        end
                    end
                    corrStruct.rois = roiNames;
                    corrStruct.corrMatrix = rCorrMatrix;
                    save(corrMatrixName{1}, 'corrStruct');
                    corrStruct.corrMatrix = zCorrMatrix;
                    save(corrMatrixName{2}, 'corrStruct'); 
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

%% Extract beta series
function meanRoi = extractBetaSeries(niiLoc, roiLoc, trimStd)
% FORMAT meanRoi = extract_beta_series(niiLoc, roiLoc, trimStd)
% Extracts and returns mean beta series from 4D niiLoc within ROI.
% Adapted from beta_series_correlation_nomars by Dennis Thompson
% without Events information.
% Calls spm_vol, roiFindIndex, adjustXyz, and spm_get_data.
%
%
% niiLoc:           Path to LSS file (4D nifti) for particular condition. 
%                   Contains data to be correlated. String.
% roiLoc:           Path to mask file to be used. String.
% trimStd:          The number of standard deviations to use if you wish 
%                   the raw data to be Windsorized. Set to 0 if the raw 
%                   data are to be used. LSS default is 3. Double.

if ~exist('trimStd', 'var'), trimStd = 3; end
threshold = 0; % Find mask values greater than this

% Get header info for beta data.
niiHeader = spm_vol(niiLoc);

% Get ROI index and transform matrix.
if exist(roiLoc, 'file')
    [XYZ, ROImat] = roiFindIndex(roiLoc, threshold);
else
    error('Mask not found: %s\n', roiLoc);
end

% Generate XYZ locations for each beta image correcting for alignment
% issues and preallocate meanRoi vector.
betaXYZ = adjustXyz(XYZ, ROImat, niiHeader);
meanRoi = zeros(1, length(betaXYZ));

% Extract mean of ROI from each beta.
for iBeta = 1:length(betaXYZ),
    betasInRoi = spm_get_data(niiHeader(iBeta), betaXYZ{iBeta});
    meanRoi(iBeta) = nanmean(betasInRoi(:));
end

clear betasInRoi

if trimStd > 0,
    meanRoi = trimTimeseries(meanRoi, trimStd);
end
end

%% Create Correlation image
function correlationMatrix = correlateBetaSeries(niiHeader, xVol, meanRoi, trimStd)
% FORMAT correlationMatrix = correlateBetaSeries(niiHeader, xVol, meanRoi, trimStd)
% Correlates meanRoi beta series with beta series of all other voxels in
% brain. Converts resulting R matrix to Z matrix.
% Adapted from beta_series_correlation_nomars by Dennis Thompson
% without Events information.
% Calls spm_vol and spm_get_data.
%
%
% niiHeader:        Header of LSS file (4D nifti) for particular condition. 
%                   Contains data to be correlated. String.
% xVol:             Select information from LSS file. Structure.
% meanRoi:          Mean beta series correlated with all other voxels' beta
%                   series in brain. Double vector.
% trimStd:          The number of standard deviations to use if you wish 
%                   the raw data to be Windsorized. Set to 0 if the raw 
%                   data are to be used. LSS default is 3. Double.

if ~exist('trimStd', 'var'), trimStd = 3; end

allDataAtBetas = spm_get_data(niiHeader, xVol.XYZ);
correlationMatrix = cell(3); correlationMatrix{1} = zeros(1, size(allDataAtBetas, 2));
for iVoxel = 1:size(allDataAtBetas, 2),
    if trimStd > 0,
        allDataAtBetas(:, iVoxel) = trimTimeseries(allDataAtBetas(:, iVoxel), trimStd);
    end
    correlationMatrix{1}(iVoxel) = corr(meanRoi, allDataAtBetas(:, iVoxel), 'type', 'Pearson');
end
correlationMatrix{2} = atanh(correlationMatrix{1});
stdZ = 1 / (length(niiHeader) - 3);
correlationMatrix{3} = correlationMatrix{2} / stdZ;
end

%% Write Correlation Image
function writeCorrelationImage(correlationMatrix, outNiftiName, xVol)
% FORMAT writeCorrelationImage(correlationMatrix, outNiftiName, xVol)
% Writes out whole brain correlation data into a nifti image. Verbatim
% Dennis Thompson's write_correlation_image.
%
%
% correlationMatrix:   A vector containing the correlations.
% outNiftiName:        What the file is to be called, can include the full path.
% xVol:                A sub set of the SPM data structure.

nvoxels = size(xVol.XYZ, 2); % length of real data in brain - should match the length of correlationMatrix
corrData = NaN * ones(xVol.DIM(1), xVol.DIM(2), xVol.DIM(3)); % array of NaN size of brain volumne

% replace NaN's with corr results
for iVoxel = 1:nvoxels,
    corrData(xVol.XYZ(1, iVoxel), xVol.XYZ(2, iVoxel), xVol.XYZ(3, iVoxel)) = correlationMatrix(iVoxel);
end

[pathstr, name, ext] = fileparts(outNiftiName);

% make nifti data structure see spm_vol - dt is [32 bit float, little endian]
corrIm = struct ('fname', [name, ext], ...
    'dim',        [xVol.DIM'], ...
    'dt',         [16, 0], ...
    'mat',        xVol.M, ...
    'pinfo',      [1 0 0]', ...
    'descript',   'beta-correlation');

cwd = pwd;
cd(pathstr);

corrIm = spm_create_vol(corrIm, 'noopen');
[~] = spm_write_vol(corrIm, corrData);

cd(cwd);
end

%% Trim Timeseries
function [y, ntrimmed] = trimTimeseries(y, sd)
% FORMAT [y, ntrimmed] = trimTimeseries(y, sd)
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
function [index, mat] = roiFindIndex(roiLoc, thresh)
% FORMAT [index, mat] = roiFindIndex(roiLoc, thresh)
% Returns the XYZ address of voxels with values greater than threshold. 
% By Dennis Thompson.
% 
%
% roiLoc:           String pointing to nifti image.
% thresh:           Threshold value, defaults to zero. Double.

if ~exist('thresh','var'),
    thresh = 0;
end

data = nifti(roiLoc);
Y = double(data.dat);
Y(isnan(Y)) = 0;
index = [];
for n = 1:size(Y,3)
    % find values greater > thresh
    [xx, yy] = find(squeeze(Y(:,:,n)) > thresh);
    if ~isempty(xx),
        zz = ones(size(xx))*n;
        index = [index [xx';yy';zz']];
    end
end

mat = data.mat;
end

%% Adjust Coordinates of ROI
function funcXyz = adjustXyz(xyz, ROImat, niiHeader)
% FORMAT funcXyz = adjustXyz(xyz, ROImat, niiHeader)
% By Dennis Thompson.
%
%
% XYZ:              Output from roiFindIndex.
% ROImat:           Output from roiFindIndex.
% niiHeader:        Header information of nifti file from spm_vol.

xyz(4,:) = 1;
funcXyz = cell(length(niiHeader));
for n = 1:length(niiHeader),
    if(iscell(niiHeader)),
       tmp = inv(niiHeader{n}.mat) * (ROImat * xyz);
    else
        tmp = inv(niiHeader(n).mat) * (ROImat * xyz);
    end
    funcXyz{n} = tmp(1:3,:);
end
end

%% Write Out Mean Image
function writeSummaryImage(cellVols, outName, expr)
% FORMAT writeSummaryImage(cellVols, outName, expr)
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
