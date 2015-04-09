function saveClustersAndEffectSize(spmFile, pThr, correction, k, maskFile)
% FORMAT saveClustersAndEffectSize(spmFile, pThr, correction, k, maskFile)
% Loops through T contrasts in SPM second-level analysis, thresholding spmT
% image based on pThr, correction, and k, and saving a mask of each cluster
% in subdirectory. Also creates a Cohen's d image for each T contrast and
% extracts mean Cohen's d for each significant cluster, summarizing in
% outputted csv.
%
%
% spmFile:          Path to SPM.mat file (including SPM.mat). String.
% pThr:             p threshold. Cell array (1x1 or 1x2) of doubles.
% correction:       Correction applied to p threshold (FWE, FDR, or unc).
%                   Cell array (same size as pThr) of strings.
% k:                Minimum acceptable cluster size. Double.
%
%
% 141019-150320 Created by Taylor Salo

%% Check inputs
if exist(spmFile, 'file')
    [path, ~] = fileparts(spmFile);
    load(spmFile);
    if isempty(path)
        path = pwd;
    end
else
    error('spmFile does not exist. Quitting.');
end

[pThr, correction, k] = checkInputs(pThr, correction, k);

%% Do everything else.
origDir = pwd;
cd(path);
fprintf('Evaluating second-level results of design: %s.\n', SPM.xsDes.Design);
if length(pThr) == 1
    fprintf(['Evaluating at p < ' num2str(pThr{1}) ' ' correction{1} ' and k > ' num2str(k) '.\n']);
    clusterExtentThresholdingDetected = false;
else
    fprintf('Evaluating at voxel-level p < %g %s and k > %d and cluster-level p < %g %s.\n', pThr{1}, correction{1}, k, pThr{2}, correction{2});
    clusterExtentThresholdingDetected = true;
    fprintf('Cluster extent thresholding detected: Only cluster-level statistics will be reported.\n');
    switch correction{2}
        case 'unc'
            clusterLevelSigCol = 4;
        case 'FWE'
            clusterLevelSigCol = 1;
        case 'FDR'
            clusterLevelSigCol = 2;
    end
end
resmsFile = SPM.VResMS.fname;
[~, ~, fileSuffix] = fileparts(resmsFile);
rpvFile = [path '/RPV' fileSuffix];

for iCon = 1:length(SPM.xCon)
    xSPM.STAT = SPM.xCon(iCon).STAT;
    xSPM.df = [SPM.xCon(iCon).eidf SPM.xX.erdf];
    xSPM.k = k;
    xSPM.VRpv = spm_vol(rpvFile);
    xSPM.M = xSPM.VRpv.mat;
    xSPM.n = 1;
    xSPM.S = SPM.xVol.S;
    xSPM.Ic = iCon;
    xSPM.R = SPM.xVol.R;
    xSPM.FWHM = SPM.xVol.FWHM;
    xSPM.DIM = SPM.xVol.DIM;
    
    spmTFile = fullfile(path, SPM.xCon(iCon).Vspm.fname);
    conName = sprintf('Contrast_%03d-%s', iCon, strrep(SPM.xCon(iCon).name, ' ', '_'));
    fprintf('\tEvaluating contrast %d, %s\n', iCon, conName);
    
    if ~isempty(maskFile)
        [~, maskName, ~] = fileparts(maskFile);
        addMask = ['_' maskName];
    else
        addMask = '';
    end
    
    if clusterExtentThresholdingDetected
        outDir = fullfile(path, sprintf('v%g_%s_k%d_c%g_%s%s_clusters/%s/',...
                                        pThr{1}, correction{1}, k, pThr{2},...
                                        correction{2}, addMask, conName));
    else
        outDir = fullfile(path, sprintf('v%g_%s_k%d%s_clusters/%s/', pThr{1},...
                                        correction{1}, k, addMask, conName));
    end
    
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    
    % Create Cohen's D image.
    VspmT = spm_vol(spmTFile);
    [contrastTValues, ~] = spm_read_vols(VspmT);
    
    if ~isempty(maskFile)
        [xyz, maskMatrix] = maskFindIndex(maskFile, 0);
        roiXyz = adjustXyz(xyz, maskMatrix, VspmT);
        ind = sub2ind(size(contrastTValues), roiXyz{1}(1, :), roiXyz{1}(2, :), roiXyz{1}(3, :))';
        maskedValues = zeros(size(contrastTValues));
        maskedValues(ind) = contrastTValues(ind);
        contrastTValues = maskedValues;
    end
    
    if strcmp(xSPM.STAT, 'T')
        dValues = (2 .* contrastTValues) ./ sqrt(xSPM.df(2));
        dFileHeader = VspmT;
        dFileHeader.fname = fullfile(outDir, ['D_' conName '.nii']);
        outDFileHeader = spm_create_vol(dFileHeader);
        spm_write_vol(outDFileHeader, dValues);
    end
    
    % Determine voxel size of T image.
    [sizeX, sizeY, sizeZ] = getVoxelSize(spmTFile);
    voxelScalar = sizeX * sizeY * sizeZ;
    xSPM.VOX = [sizeX sizeY sizeZ];
    xSPM.units = {'mm' 'mm' 'mm'};
    
    % Set up output csv.
    clear outStruct
    outStruct{1}.header{1} = 'Region of Activation'; outStruct{2}.header{1} = 'BA';
    outStruct{3}.header{1} = 'L/R'; outStruct{4}.header{1} = 'k (mm3)';
    outStruct{7}.header{1} = 'x'; outStruct{8}.header{1} = 'y'; outStruct{9}.header{1} = 'z';
    
    if clusterExtentThresholdingDetected
        outStruct{5}.header{1} = 'Mean T';
        outStruct{6}.header{1} = 'Mean D';
    else
        outStruct{5}.header{1} = 'T';
        outStruct{6}.header{1} = 'D';
    end
    
    roaCol = 1; baCol = 2; lrCol = 3; kCol = 4; tCol = 5; dCol = 6; xCol = 7; yCol = 8; zCol = 9;
    outStruct{1}.col{1} = ''; outStruct{2}.col{1} = ''; outStruct{3}.col{1} = '';
    outStruct{4}.col{1} = ''; outStruct{5}.col{1} = ''; outStruct{6}.col{1} = '';
    outStruct{7}.col{1} = ''; outStruct{8}.col{1} = ''; outStruct{9}.col{1} = '';
    
    % Create masks of all significant clusters and determine mean
    % Cohen's D of each cluster.
    switch correction{1}
        case 'unc'
            xSPM.u = spm_u(pThr{1}, xSPM.df, xSPM.STAT);
        case 'FWE'
            xSPM.u  = spm_uc(pThr{1}, xSPM.df, xSPM.STAT, xSPM.R, 1, xSPM.S);
        case 'FDR'
            xSPM.u  = spm_uc_FDR(pThr{1}, xSPM.df, xSPM.STAT, 1, VspmT, 0);
        otherwise
            error('Variable corr must be either ''unc'', ''FWE'', or ''FDR''.');
    end
    
    clusterFileHeader = VspmT;
    allClustersFileHeader = VspmT;
    [allClustVals, nClusters] = spm_bwlabel(double(contrastTValues > xSPM.u), 18);
    [clustSize, clustNum] = sort(histc(allClustVals(:), 0:nClusters), 1, 'descend');
    clustSize = clustSize(2:end); clustNum = clustNum(2:end) - 1;
    clustNum = clustNum(clustSize >= k); clustSize = clustSize(clustSize >= k);
    
    [x, y, z] = ind2sub(size(contrastTValues), find(contrastTValues > xSPM.u));
    XYZ = [x'; y'; z'];
    Z = min(Inf, spm_get_data(SPM.xCon(iCon).Vspm, XYZ));
    
    V2R = 1/prod(SPM.xVol.FWHM(SPM.xVol.DIM > 1));
    
    if strcmp(xSPM.STAT, 'T') || ~isempty(strfind(which('spm'), 'spm12'))
        continueClusterExtentThresholding = true;
        [uc, xSPM.Pc, ue] = spm_uc_clusterFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, Z, XYZ, V2R, xSPM.u);
    else
        continueClusterExtentThresholding = false;
        uc = NaN; ue = NaN; xSPM.Pc = ones(1, nClusters) * .05;
    end
    
    if clusterExtentThresholdingDetected && ~continueClusterExtentThresholding
        fprintf(['\tAttempting to use cluster-extent thresholding with F-contrast.\n'...
            '\t\tThis is possible in SPM12, but not the currently loaded version of SPM.\n']);
        continue
    end
    
    [up, xSPM.Pp] = spm_uc_peakFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, Z, XYZ, xSPM.u);
    
    A = spm_clusters(XYZ);
    Q = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= xSPM.k
            Q = [Q j];
        end
    end
    
    xSPM.Z = Z(:, Q);
    xSPM.XYZ = XYZ(:, Q);
    xSPM.XYZmm = xSPM.XYZ;
    
    for jVox = 1:size(xSPM.XYZ, 2)
        xSPM.XYZmm(:, jVox) = (xSPM.XYZ(:, jVox).' * xSPM.M(1:3, 1:3) + xSPM.M(1:3, 4)').';
    end
    uu = spm_uc(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.S);
    
    xSPM.uc = [uu up ue uc];
    
    table = spm_list_edited(xSPM, 5, 8);
    
    % Cluster row index.
    clustLevelCol = table(:, 3);
    clustIdx = find(~cellfun(@isempty, clustLevelCol));
    
    % Pare down table to cluster-level significant clusters if cluster-
    % extent based thresholding is being applied. Else, keep as is.
    if clusterExtentThresholdingDetected
        for iClust = length(clustIdx):-1:1
            if iClust == length(clustIdx)
                if table{clustIdx(iClust), clusterLevelSigCol} > pThr{2}
                    table(clustIdx(iClust):end, :) = [];
                end
            else
                if table{clustIdx(iClust), clusterLevelSigCol} > pThr{2}
                    table(clustIdx(iClust):clustIdx(iClust + 1) - 1, :) = [];
                end
            end
        end
    end
    
    clustLevelCol = table(:, 3);
    clustIdx = find(~cellfun(@isempty, clustLevelCol));
    clear x y z XYZ
    
    for jClust = 1:length(clustNum)
        oneClustVals = zeros(size(allClustVals));
        oneClustVals(allClustVals == clustNum(jClust)) = 1;
        [x, y, z] = ind2sub(size(oneClustVals), find(oneClustVals == 1));
        clusterCoordinates = [x'; y'; z'];
        clustNumber = 0;
        
        % Determine which cluster you're looking at and fill in that
        % position in the output table. Basically we're merging the
        % information about the clusters from the table with the
        % information found through more direct means, including effect
        % size.
        for kClust = 1:length(clustIdx)
            peakMM(1:3) = (table{clustIdx(kClust), 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);
            if any(ismember(clusterCoordinates.', peakMM, 'rows'))
                clustNumber = kClust;
                clustPeakMM = peakMM;
            end
        end
        
        % Skip clusters not found in the table (because they're not
        % significant). Only skips clusters if cluster-extent
        % thresholding was used.
        if clustNumber ~= 0
            if ~exist('allClustMask', 'var')
                allClustMask = oneClustVals;
            else
                allClustMask = allClustMask + oneClustVals;
            end
            peakCoord = table{clustIdx(clustNumber), 10}.';
            clusterFile = fullfile(outDir, sprintf('Cluster_%03d_%d_%d_%d.nii', clustNumber, peakCoord(1), peakCoord(2), peakCoord(3)));
            clusterFileHeader.fname = clusterFile;
            spm_write_vol(clusterFileHeader, oneClustVals);
            
            % Fill in output csv.
            if peakCoord(1) < 0
                outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'L';
            elseif peakCoord(1) > 0
                outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'R';
            else
                outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'I';
            end
            
            outStruct{kCol}.col{clustIdx(clustNumber), 1} = clustSize(jClust) * voxelScalar;
            
            outStruct{xCol}.col{clustIdx(clustNumber), 1} = peakCoord(1);
            outStruct{yCol}.col{clustIdx(clustNumber), 1} = peakCoord(2);
            outStruct{zCol}.col{clustIdx(clustNumber), 1} = peakCoord(3);
            
            if clusterExtentThresholdingDetected
                if strcmp(xSPM.STAT, 'T')
                    outStruct{dCol}.col{clustIdx(clustNumber), 1} = mean(dValues(find(oneClustVals == 1)));
                else
                    outStruct{dCol}.col{clustIdx(clustNumber), 1} = '';
                end
                [outString, baString] = getBaAndRegion(clusterFile);
                outStruct{roaCol}.col{clustIdx(clustNumber), 1} = outString;
                outStruct{baCol}.col{clustIdx(clustNumber), 1} = baString;
                outStruct{tCol}.col{clustIdx(clustNumber), 1} = mean(contrastTValues(find(oneClustVals == 1)));
            else
                if strcmp(xSPM.STAT, 'T')
                    outStruct{dCol}.col{clustIdx(clustNumber), 1} = dValues(clustPeakMM(1), clustPeakMM(2), clustPeakMM(3));
                else
                    outStruct{dCol}.col{clustIdx(clustNumber), 1} = '';
                end
                [outString, baString] = getBaAndRegion(peakCoord);
                outStruct{roaCol}.col{clustIdx(clustNumber), 1} = outString;
                outStruct{baCol}.col{clustIdx(clustNumber), 1} = baString;
                outStruct{tCol}.col{clustIdx(clustNumber), 1} = num2str(table{clustIdx(clustNumber), 7});
            end
            
            if clustNumber == length(clustIdx)
                [nRows, ~] = size(table);
            else
                nRows = clustIdx(clustNumber + 1) - 1;
            end
            
            for mSubClust = clustIdx(clustNumber) + 1:nRows
                if table{mSubClust, 10}(1) < 0
                    outStruct{lrCol}.col{mSubClust, 1} = 'L';
                elseif table{mSubClust, 10}(1) > 0
                    outStruct{lrCol}.col{mSubClust, 1} = 'R';
                else
                    outStruct{lrCol}.col{mSubClust, 1} = 'I';
                end
                
                outStruct{kCol}.col{mSubClust, 1} = '';
                subPeakMm = (table{mSubClust, 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);
                
                if clusterExtentThresholdingDetected
                    outStruct{dCol}.col{mSubClust, 1} = '';
                    outStruct{tCol}.col{mSubClust, 1} = '';
                    outStruct{roaCol}.col{mSubClust, 1} = '';
                    outStruct{baCol}.col{mSubClust, 1} = '';
                else
                    if strcmp(xSPM.STAT, 'T')
                        outStruct{dCol}.col{mSubClust, 1} = dValues(subPeakMm(1), subPeakMm(2), subPeakMm(3));
                    else
                        outStruct{dCol}.col{mSubClust, 1} = '';
                    end
                    [outString, baString] = getBaAndRegion(table{mSubClust, 10}');
                    outStruct{roaCol}.col{mSubClust, 1} = ['\t' outString];
                    outStruct{baCol}.col{mSubClust, 1} = baString;
                    outStruct{tCol}.col{mSubClust, 1} = num2str(table{mSubClust, 7});
                end
                
                outStruct{xCol}.col{mSubClust, 1} = table{mSubClust, 10}(1);
                outStruct{yCol}.col{mSubClust, 1} = table{mSubClust, 10}(2);
                outStruct{zCol}.col{mSubClust, 1} = table{mSubClust, 10}(3);
            end
        end
    end
    if exist('allClustMask', 'var')
        allSignificantClusterValues = contrastTValues .* allClustMask;
        allClustersFileHeader.fname = fullfile(outDir, 'allClusterMask.nii');
        spm_write_vol(allClustersFileHeader, allClustMask);
        
        allClustersFileHeader.fname = fullfile(outDir, 'allClusterVals.nii');
        spm_write_vol(allClustersFileHeader, allSignificantClusterValues);
        clear allClustMask
    end
    
    fprintf('\t\t%d out of %d clusters are larger than %d voxels.\n', length(clustIdx), nClusters, k);
    writeCsv(outStruct, fullfile(outDir, sprintf('ClusterReport_%03d.csv', length(clustIdx))));
end
cd(origDir);
fprintf('Done.\n\n');
end

%% Check Inputs
function [pThr, correction, k] = checkInputs(pThr, correction, k)
% FORMAT [pThr, correction, k] = checkInputs(pThr, correction, k)
% By Taylor Salo.
if iscell(pThr)
    [m, n] = size(pThr);
    if m * n > 2
        fprintf('pThr is %dx%d.\n\tSetting pThr to {0.01 0.05} and corr to {''unc'' ''FWE''}.\n', m, n);
        pThr = {0.01 0.05};
        correction = {'unc' 'FWE'};
    elseif m * n == 2
        for iP = 1:length(pThr)
            if ~isa(pThr{iP}, 'double')
                fprintf('pThr{%d} is not double.\n\tSetting pThr to {0.01 0.05} and corr to {''unc'' ''FWE''}.\n', iP);
                pThr = {0.01 0.05};
                correction = {'unc' 'FWE'};
                break
            end
        end
    else
        if ~isa(pThr{1}, 'double')
            fprintf('pThr{1} is not double.\n\tSetting pThr to {0.001} and corr to {''unc''}.\n');
            pThr = {0.001};
            correction = {'unc'};
        end
    end
else
    fprintf('pThr is not a 1x2 or 2x1 cell array of doubles.\n\tSetting pThr to {0.01 0.05} and corr to {''unc'' ''FWE''}.\n');
    pThr = {0.01 0.05};
    correction = {'unc' 'FWE'};
end
if iscellstr(correction)
    if length(correction) ~= length(pThr)
        fprintf('pThr and corr are different lengths.\n');
        if length(pThr) == 2
            fprintf('\tSetting corr to {''unc'' ''FWE''}.\n');
            correction = {'unc' 'FWE'};
        else
            fprintf('\tSetting corr to {''unc''}.\n');
            correction = {'unc'};
        end
    else
        if length(correction) == 1
            if ~cellstrfind(correction{1}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of ''unc'', ''FWE'', or ''FDR''.\n\tSetting corr to {''unc''}.\n');
                correction = {'unc'};
            end
        else
            if ~cellstrfind(correction{1}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of ''unc'', ''FWE'', or ''FDR''.\n\tSetting corr{1} to ''unc''.\n');
                correction{1} = 'unc';
            end
            if ~cellstrfind(correction{2}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of ''unc'', ''FWE'', or ''FDR''.\n\tSetting corr{2} to ''FWE''.\n');
                correction{2} = 'FWE';
            end
        end
    end
else
    fprintf('corr is not a cellstr.\n');
    if length(pThr) == 2
        fprintf('\tSetting corr to {''unc'' ''FWE''}.\n');
        correction = {'unc' 'FWE'};
    else
        fprintf('\tSetting corr to {''unc''}.\n');
        correction = {'unc'};
    end
end
if ~isa(k, 'double')
    fprintf('Warning, your set k is not a double. Setting to default 5.\n');
    k = 5;
end
end

%% Determine Voxel Size in Nifti File
function [sizeX, sizeY, sizeZ] = getVoxelSize(niiFile)
% FORMAT [sizeX, sizeY, sizeZ] = getVoxelSize(niiFile)
% Determines size of voxel in X, Y, and Z dimensions. This assumes that the
% matrix is not diagonal, but works fine if it is.
% By Taylor Salo.
%
%
% niiFile:          Nifti file from which voxel size will be determined.
%                   String.
%
% sizeX (output):   Size in mm of voxel in X dimension. Double.
% sizeY (output):   Size in mm of voxel in Y dimension. Double.
% sizeZ (output):   Size in mm of voxel in Z dimension. Double.
V = spm_vol(niiFile);
mniO = (V.mat * [1 1 1 1]')';
mniX = (V.mat * [2 1 1 1]')';
mniY = (V.mat * [1 2 1 1]')';
mniZ = (V.mat * [1 1 2 1]')';
sizeX = pdist([mniO; mniX], 'euclidean');
sizeY = pdist([mniO; mniY], 'euclidean');
sizeZ = pdist([mniO; mniZ], 'euclidean');
end

%% Extract Coordinates of ROI
function [index, mat] = maskFindIndex(roiLoc, thresh)
% FORMAT [index, mat] = maskFindIndex(roiLoc, thresh)
% Returns the XYZ address of voxels with values greater than threshold.
% By Dennis Thompson.
%
% Inputs:
% roiLoc:           String pointing to nifti image.
% thresh:           Threshold value, defaults to zero. Double.
%
% Outputs:
% index:            Index of ROI XYZ coordinates.
% mat:              Mask affine transformation matrix.
%
%
% 080101 Created by Dennis Thompson

if ~exist('thresh','var'),
    thresh = 0;
end

data = nifti(roiLoc);
Y = double(data.dat);
Y(isnan(Y)) = 0;
index = [];
for iSheet = 1:size(Y, 3)
    % find values greater > thresh
    [xx, yy] = find(squeeze(Y(:, :, iSheet)) > thresh);
    if ~isempty(xx)
        zz = ones(size(xx)) * iSheet;
        index = [index [xx'; yy'; zz']];
    end
end

mat = data.mat;
end

%% Adjust Coordinates of ROI
function funcXYZ = adjustXyz(XYZ, ROImat, V)
% FORMAT funcXYZ = adjustXyz(XYZ, ROImat, V)
% By Dennis Thompson.
%
%
% XYZ:              XYZ coordinates of ones in binary matrix.
% ROImat:           Transformation matrix from ROI file.
% V:                Header information of nifti file from spm_vol.
XYZ(4, :) = 1;
funcXYZ = cell(length(V));
for n = 1:length(V)
    if iscell(V)
        tmp = inv(V{n}.mat) * (ROImat * XYZ);
    else
        tmp = inv(V(n).mat) * (ROImat * XYZ);
    end
    funcXYZ{n} = tmp(1:3, :);
end
end
