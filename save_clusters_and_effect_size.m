function save_clusters_and_effect_size(spmFile, pThr, corr, k)
% FORMAT save_clusters_and_effect_size(spmFile, pThr, corr, k)
% Loops through T contrasts in SPM second-level analysis, thresholding spmT
% image based on pThr, corr, and k, and saving a mask of each cluster in
% subdirectory. Also creates a Cohen's d image for each T contrast and
% extracts mean Cohen's d for each significant cluster, summarizing in
% outputted csv.
%
%
% spmFile: Path to SPM.mat file (including SPM.mat). String.
% pThr:    p threshold. Double.
% corr:    Correction applied to p threshold (FWE, FDR, or unc). String.
% k:       Minimum acceptable cluster size. Double.

if exist(spmFile, 'file')
    [path, ~] = fileparts(spmFile);
    load(spmFile);
else
    error('spmFile does not exist. Quitting.');
end

if ischar(pThr)
    fprintf('Warning, your set p threshold is a string. Trying to convert to double.\n')
    pThr = str2double(pThr);
    if isnan(pThr)
        fprintf('P threshold could not be converted to double, setting to default 0.001.\n');
        pThr = 0.001;
    end
elseif ~isa(pThr, 'double')
    fprintf('Warning, your set p threshold is not a double or a string. Setting to default 0.001.\n');
    pThr = 0.001;
end

if ischar(k)
    fprintf('Warning, your set k is a string. Trying to convert to double.\n')
    k = str2double(k);
    if isnan(k)
        fprintf('P threshold could not be converted to double, setting to default 5.\n');
        k = 5;
    end
elseif ~isa(k, 'double')
    fprintf('Warning, your set k is not a double or a string. Setting to default 5.\n');
    k = 5;
end

design = SPM.xsDes.Design;
fprintf(['Evaluating second-level results of design: ' design '.\n']);
fprintf(['Evaluating at p < ' num2str(pThr) ' ' corr ' and k > ' num2str(k) '.\n']);

for iCon = 1:length(SPM.xCon)
    STAT = SPM.xCon(iCon).STAT;

    if strcmp(STAT, 'T')
        spmT = [path '/' SPM.xCon(iCon).Vspm.fname];
        conName = ['Contrast_' sprintf('%03d', iCon) '-' strrep(SPM.xCon(iCon).name, ' ', '_')];
        fprintf(['\tEvaluating contrast ' num2str(iCon) ', ' conName '\n']);
        outDir = [path '/' num2str(pThr) '_' corr '_clusters/' conName '/'];

        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        % Create Cohen's D image.
        VspmT = spm_vol(spmT);
        [Tvals, ~] = spm_read_vols(VspmT);
        df = [SPM.xCon(iCon).eidf SPM.xX.erdf];
        Dvals = (2 .* Tvals) ./ sqrt(df(2));
        dHeader = VspmT;
        dHeader.fname = [outDir 'D_' conName '.nii'];
        outDHeader = spm_create_vol(dHeader);
        spm_write_vol(outDHeader, Dvals(:, :, :));

        % Determine voxel size of T image.
        [sizeX, sizeY, sizeZ] = get_voxel_size(spmT);
        voxelScalar = sizeX * sizeY * sizeZ;

        % Set up output csv.
        outStruct{1}.header{1} = 'x'; outStruct{2}.header{1} = 'y'; outStruct{3}.header{1} = 'z';
        outStruct{4}.header{1} = 'k'; outStruct{5}.header{1} = 'mm3'; outStruct{6}.header{1} = 'peak_T';
        outStruct{7}.header{1} = 'peak_D'; outStruct{8}.header{1} = 'mean_D';
        outStruct{1}.col{1} = ''; outStruct{2}.col{1} = ''; outStruct{3}.col{1} = '';
        outStruct{4}.col{1} = ''; outStruct{5}.col{1} = ''; outStruct{6}.col{1} = '';
        outStruct{7}.col{1} = ''; outStruct{8}.col{1} = '';

        % Create masks of all significant clusters and determine mean
        % Cohen's D of each cluster.
        switch corr
            case 'unc'
                tThr = spm_u(pThr, df, STAT);
            case 'FWE'
                tThr  = spm_uc(pThr, df, STAT, SPM.xVol.R, 1, SPM.xVol.S);
            case 'FDR'
                tThr  = spm_uc_FDR(pThr, df, STAT, 1, VspmT, 0);
            otherwise
                error('Variable corr must be either "unc", "FWE", or "FDR".');
        end

        clustHeader = VspmT;
        rawData = spm_read_vols(VspmT);
        [allClustMat, nClusters] = spm_bwlabel(double(rawData > tThr), 18);

        [clustSize, clustNum] = sort(histc(allClustMat(:), 0:nClusters), 1, 'descend');
        clustSize = clustSize(2:end); clustNum = clustNum(2:end) - 1;
        clustNum = clustNum(clustSize >= k); clustSize = clustSize(clustSize >= k);

        for jClust = 1:length(clustSize)
            oneClustMat = zeros(size(allClustMat));
            oneClustMat(allClustMat == clustNum(jClust)) = 1;
            [x, y, z] = ind2sub(size(oneClustMat), find(oneClustMat == 1));
            XYZ = [x'; y'; z'];
            betaXYZ = adjust_XYZ(XYZ, clustHeader.mat, dHeader);
            roi = spm_get_data(dHeader.fname, betaXYZ{1});
            roi(isnan(roi)) = 0;
            [peakValue, peakCoord, peakMM] = get_peak_d(VspmT.fname, betaXYZ);
            clustHeader.fname = [outDir 'Cluster_' sprintf('%03d', jClust) '_' num2str(peakCoord(1)) '_' num2str(peakCoord(2)) '_' num2str(peakCoord(3)) '.nii'];
            spm_write_vol(clustHeader, oneClustMat);

            % Fill in output csv.
            outStruct{1}.col{jClust, 1} = peakCoord(1);
            outStruct{2}.col{jClust, 1} = peakCoord(2);
            outStruct{3}.col{jClust, 1} = peakCoord(3);
            outStruct{4}.col{jClust, 1} = clustSize(jClust);
            outStruct{5}.col{jClust, 1} = clustSize(jClust) * voxelScalar;
            outStruct{6}.col{jClust, 1} = peakValue;
            outStruct{7}.col{jClust, 1} = Dvals(peakMM(1), peakMM(2), peakMM(3));
            outStruct{8}.col{jClust, 1} = sum(roi(:)) / length(find(roi));
        end

        fprintf('\t\t%d out of %d clusters are larger than %d voxels.\n', length(clustSize), nClusters, k);
        write_csv(outStruct, [outDir 'ClusterReport_' sprintf('%03d', length(clustSize)) '.csv']);
    else
        fprintf(['\tSkipping contrast ' num2str(iCon) ', ' SPM.xCon(iCon).name ', because it is an F con.\n']);
    end
end

fprintf('Done.\n\n');
end

%% Determine Voxel Size in Nifti File
function [sizeX, sizeY, sizeZ] = get_voxel_size(niiFile)
% FORMAT [sizeX, sizeY, sizeZ] = get_voxel_size(niiFile)
% Determines size of voxel in X, Y, and Z dimensions. This assumes that the
% matrix is not diagonal, but works fine if it is.
V = spm_vol(niiFile);
mniO = (V.mat * [1 1 1 1]')';
mniX = (V.mat * [2 1 1 1]')';
mniY = (V.mat * [1 2 1 1]')';
mniZ = (V.mat * [1 1 2 1]')';
sizeX = pdist([mniO; mniX], 'euclidean');
sizeY = pdist([mniO; mniY], 'euclidean');
sizeZ = pdist([mniO; mniZ], 'euclidean');
end

%% Extract Value, MNI Coordinates, and Voxel Coordinates of Peak Value in ROI.
function [peakValue, peakCoord, peakMM] = get_peak_d(niiFile, betaXYZ)
% FORMAT [peakValue, peakCoord, peakMM] = get_peak_d(niiFile, betaXYZ)
% Extracts coordinates, value, and location in mm of peak value in roi from
% nifti image.
% Extract values from ROI and get max and location in XYZ and coordinates.
V = spm_vol(niiFile);
valuesInROI = spm_get_data(V(1), betaXYZ{1});
[peakValue, idx] = max(valuesInROI);
maxLoc = betaXYZ{1}(:, idx)';
peakCoord(1:3) = maxLoc * V.mat(1:3, 1:3) + V.mat(1:3, 4)';
peakMM = uint8(maxLoc);
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
