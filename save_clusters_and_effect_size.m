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

[path, ~] = fileparts(spmFile);
load(spmFile);

design = SPM.xsDes.Design;

conCounter = 1;
maxClusters = 1;

fprintf(['Evaluating second-level results of design: ' design '.\n']);
fprintf(['Evaluating at p < ' num2str(pThr) ' ' corr ' and k > ' num2str(k) '.\n']);

for iCon = 1:length(SPM.xCon)
    STAT = SPM.xCon(iCon).STAT;
    
    if strcmp(STAT, 'T')
        spmT = [path '/' SPM.xCon(iCon).Vspm.fname];
        
        conName = SPM.xCon(iCon).name;
        conName = strrep(conName, ' ', '_');
        fprintf(['\tEvaluating contrast ' num2str(iCon) ', ' conName '\n']);
        outDir = [path '/' num2str(pThr) '_' corr '_clusters/' conName '/'];

        % Create Cohen's D image.
        [VspmT, Tvals] = load_nii_spm(spmT);
        df = [SPM.xCon(iCon).eidf SPM.xX.erdf];
        D_name = 'D_';
        Dvals = (2 .* Tvals) ./ sqrt(df(2));
        dHeader = VspmT;
        dHeader.fname = [outDir D_name conName '.nii'];
        
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        
        save_nii_spm(dHeader, Dvals);
        
        outStruct{conCounter}.header{1} = [conName '_X'];
        outStruct{conCounter+1}.header{1} = [conName '_Y'];
        outStruct{conCounter+2}.header{1} = [conName '_Z'];
        outStruct{conCounter+3}.header{1} = [conName '_size'];
        outStruct{conCounter+4}.header{1} = [conName '_peakT'];
        outStruct{conCounter+5}.header{1} = [conName '_peakD'];
        outStruct{conCounter+6}.header{1} = [conName '_meanD'];
        outStruct{conCounter}.col{1} = '';
        outStruct{conCounter+1}.col{1} = '';
        outStruct{conCounter+2}.col{1} = '';
        outStruct{conCounter+3}.col{1} = '';
        outStruct{conCounter+4}.col{1} = '';
        outStruct{conCounter+5}.col{1} = '';
        outStruct{conCounter+6}.col{1} = '';
        
        % Create masks of all clusters and determine mean Cohen's D of each cluster.
        switch corr
            case 'unc'
                tThr = spm_u(pThr, df, STAT);
            case 'FWE'
                tThr  = spm_uc(pThr, df, STAT, SPM.xVol.R, 1, SPM.xVol.S);
            case 'FDR'
                tThr  = spm_uc_FDR(pThr, df, STAT, 1, VspmT, 0);
        end
        
        outHeader = VspmT;
        rawData = spm_read_vols(VspmT);
        [allClustMat, nClusters] = spm_bwlabel(double(rawData > tThr), 18);

        [clustSize, clustNum] = sort(histc(allClustMat(:), 0:nClusters), 1, 'descend');
        clustSize = clustSize(2:end); clustNum = clustNum(2:end) - 1;
        clustNum = clustNum(clustSize >= k); clustSize = clustSize(clustSize >= k);
        
        for jClust = 1:length(clustSize)
            oneClustMat = zeros(size(allClustMat));
            oneClustMat(allClustMat == clustNum(jClust)) = 1;
            outHeader.fname = [outDir 'Cluster_' sprintf('%03d', jClust) '.nii'];
            spm_write_vol(outHeader, oneClustMat);
            
            [XYZ ROImat] = roi_find_index(outHeader.fname, 0);
            betaXYZ = adjust_XYZ(XYZ, ROImat, dHeader);
            roi = spm_get_data(dHeader.fname, betaXYZ{1});
            roi(isnan(roi)) = 0;
            
            [peakCoord, peakValue, peakMM] = get_peak_d(dHeader.fname, outHeader.fname);
            
            outStruct{conCounter}.col{jClust, 1} = peakCoord(1);
            outStruct{conCounter+1}.col{jClust, 1} = peakCoord(2);
            outStruct{conCounter+2}.col{jClust, 1} = peakCoord(3);
            outStruct{conCounter+3}.col{jClust, 1} = clustSize(jClust);
            outStruct{conCounter+4}.col{jClust, 1} = Tvals(peakMM(1), peakMM(2), peakMM(3));
            outStruct{conCounter+5}.col{jClust, 1} = peakValue;
            outStruct{conCounter+6}.col{jClust, 1} = sum(roi(:)) / length(find(roi));
        end
        
        if length(clustSize) > maxClusters
            maxClusters = length(clustSize);
        end
        
        fprintf('\t\t%d out of %d clusters are larger than %d voxels.\n', length(clustSize), nClusters, k);
        conCounter = conCounter + 7;
    else
        fprintf(['\tSkipping contrast ' num2str(iCon) ', ' SPM.xCon(iCon).name ', because it is an F con.\n']);
    end
end

% Write out summary csv.
outMainDir = fileparts(outDir(1:end-1));
if length(outStruct{1}.col) < maxClusters
    outStruct{1}.col{maxClusters} = '';
end

write_struct(outStruct, [outMainDir '/clusters.csv']);
fprintf('Done.\n\n');
end

%%
function [peakCoord, peakValue, peakMM] = get_peak_d(niiFile, roiFile)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

V = spm_vol(niiFile);

% Get ROI index and transform matrix.
if exist(roiFile, 'file')
    [XYZ ROImat] = roi_find_index(roiFile);
else
    error('Mask not found: %s\n', roiFile);
end

% Generate XYZ locations for each beta image correcting for alignment
% issues and preallocate mean_roi vector.
betaXYZ = adjust_XYZ(XYZ, ROImat, V);

% Extract mean of ROI from each beta.
betasInROI = spm_get_data(V(1), betaXYZ{1});
[peakValue, idx] = max(betasInROI);
maxLoc = betaXYZ{1}(:, idx);
maxLocReal = [26 -60 -24];

peakCoord = mm2vox(maxLoc, niiFile);
peakMM = maxLoc.';
end

%%
function vox = mm2vox(maxLoc, niiFile)
v2m = spm_get_space(niiFile);
m2v = inv(v2m);

vox(1:3) = maxLoc' * v2m(1:3, 1:3) + v2m(1:3, 4)';
end
