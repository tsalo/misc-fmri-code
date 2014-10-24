function saveClustersAndEffectSize(spmFile, pThr, corr, k)
% FORMAT saveClustersAndEffectSize(spmFile, pThr, corr, k)
% Loops through T contrasts in SPM second-level analysis, thresholding spmT
% image based on pThr, corr, and k, and saving a mask of each cluster in
% subdirectory. Also creates a Cohen's d image for each T contrast and
% extracts mean Cohen's d for each significant cluster, summarizing in
% outputted csv.
%
%
% spmFile:          Path to SPM.mat file (including SPM.mat). String.
% pThr:             p threshold. Cell array (1x1 or 1x2) of doubles.
% corr:             Correction applied to p threshold (FWE, FDR, or unc).
%                   Cell array (same size as pThr) of strings.
% k:                Minimum acceptable cluster size. Double.

%% Check inputs
if exist(spmFile, 'file')
    [path, ~] = fileparts(spmFile);
    load(spmFile);
else
    error('spmFile does not exist. Quitting.');
end
if iscell(pThr)
    [m, n] = size(pThr);
    if m * n > 2
        fprintf(['pThr is ' num2str(m) 'x' num2str(n) '.\n\tSetting pThr to {0.01 0.05} and corr to {"unc" "FWE"}.\n']);
        pThr = {0.01 0.05};
        corr = {'unc' 'FWE'};
    elseif m * n == 2
        for iP = 1:length(pThr)
            if ~isa(pThr{iP}, 'double')
                fprintf(['pThr{' num2str(iP) '} is not double.\n\tSetting pThr to {0.01 0.05} and corr to {"unc" "FWE"}.\n']);
                pThr = {0.01 0.05};
                corr = {'unc' 'FWE'};
                break
            end
        end
    else
        if ~isa(pThr{1}, 'double')
            fprintf('pThr{1} is not double.\n\tSetting pThr to {0.001} and corr to {"unc"}.\n');
            pThr = {0.001};
            corr = {'unc'};
        end
    end
else
    fprintf('pThr is not a 1x2 or 2x1 cell array of doubles.\n\tSetting pThr to {0.01 0.05} and corr to {"unc" "FWE"}.\n');
    pThr = {0.01 0.05};
    corr = {'unc' 'FWE'};
end
if iscellstr(corr)
    if length(corr) ~= length(pThr)
        fprintf('pThr and corr are different lengths.\n');
        if length(pThr) == 2
            fprintf('\tSetting corr to {"unc" "FWE"}.\n');
            corr = {'unc' 'FWE'};
        else
            fprintf('\tSetting corr to {"unc"}.\n');
            corr = {'unc'};
        end
    else
        if length(corr) == 1
            if ~cellstrfind(corr{1}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of "unc", "FWE" or, "FDR".\n\tSetting corr to {"unc"}.\n');
                corr = {'unc'};
            end
        else
            if ~cellstrfind(corr{1}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of "unc", "FWE" or, "FDR".\n\tSetting corr{1} to "unc".\n');
                corr{1} = 'unc';
            end
            if ~cellstrfind(corr{2}, {'unc' 'FWE' 'FDR'}, 'exact')
                fprintf('corr must be composed of "unc", "FWE" or, "FDR".\n\tSetting corr{2} to "FWE".\n');
                corr{2} = 'FWE';
            end
        end
    end
else
    fprintf('corr is not a cellstr.\n');
    if length(pThr) == 2
        fprintf('\tSetting corr to {"unc" "FWE"}.\n');
        corr = {'unc' 'FWE'};
    else
        fprintf('\tSetting corr to {"unc"}.\n');
        corr = {'unc'};
    end
end
if ~isa(k, 'double')
    fprintf('Warning, your set k is not a double. Setting to default 5.\n');
    k = 5;
end

%% Do everything else.
origDir = pwd;
cd(path);
design = SPM.xsDes.Design;
fprintf(['Evaluating second-level results of design: ' design '.\n']);
if length(pThr) == 1
    fprintf(['Evaluating at p < ' num2str(pThr{1}) ' ' corr{1} ' and k > ' num2str(k) '.\n']);
else
    fprintf(['Evaluating at voxel-level p < ' num2str(pThr{1}) ' ' corr{1} ' and k > ' num2str(k) ' and cluster-level p < ' num2str(pThr{2}) ' ' corr{2} '.\n']);
end
rpvFile = [path '/RPV.img'];

for iCon = 1:length(SPM.xCon)
    xSPM.STAT = SPM.xCon(iCon).STAT;

    if strcmp(xSPM.STAT, 'T')
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
        
        spmT = [path '/' SPM.xCon(iCon).Vspm.fname];
        conName = ['Contrast_' sprintf('%03d', iCon) '-' strrep(SPM.xCon(iCon).name, ' ', '_')];
        fprintf(['\tEvaluating contrast ' num2str(iCon) ', ' conName '\n']);
        
        if length(pThr) == 1
            outDir = [path '/v' num2str(pThr{1}) '_' corr{1} '_k' num2str(k) '_clusters/' conName '/'];
        else
            outDir = [path '/v' num2str(pThr{1}) '_' corr{1} '_k' num2str(k) '_c' num2str(pThr{2}) '_' corr{2} '_clusters/' conName '/'];
        end
        
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        % Create Cohen's D image.
        VspmT = spm_vol(spmT);
        [Tvals, ~] = spm_read_vols(VspmT);
        Dvals = (2 .* Tvals) ./ sqrt(xSPM.df(2));
        dHeader = VspmT;
        dHeader.fname = [outDir 'D_' conName '.nii'];
        outDHeader = spm_create_vol(dHeader);
        spm_write_vol(outDHeader, Dvals(:, :, :));

        % Determine voxel size of T image.
        [sizeX, sizeY, sizeZ] = get_voxel_size(spmT);
        voxelScalar = sizeX * sizeY * sizeZ;
        xSPM.VOX = [sizeX sizeY sizeZ];
        xSPM.units = {'mm' 'mm' 'mm'};
        
        % Set up output csv.
        clear outStruct
        outStruct{1}.header{1} = 'Region of Activation'; outStruct{2}.header{1} = 'BA'; outStruct{3}.header{1} = 'L/R';
        outStruct{4}.header{1} = 'k (mm3)'; outStruct{5}.header{1} = 'T'; outStruct{6}.header{1} = 'D';
        outStruct{7}.header{1} = 'x'; outStruct{8}.header{1} = 'y'; outStruct{9}.header{1} = 'z';
        roaCol = 1; baCol = 2; lrCol = 3; kCol = 4; tCol = 5; dCol = 6; xCol = 7; yCol = 8; zCol = 9;
        outStruct{1}.col{1} = ''; outStruct{2}.col{1} = ''; outStruct{3}.col{1} = '';
        outStruct{4}.col{1} = ''; outStruct{5}.col{1} = ''; outStruct{6}.col{1} = '';
        outStruct{7}.col{1} = ''; outStruct{8}.col{1} = ''; outStruct{9}.col{1} = '';

        % Create masks of all significant clusters and determine mean
        % Cohen's D of each cluster.
        switch corr{1}
            case 'unc'
                xSPM.u = spm_u(pThr{1}, xSPM.df, xSPM.STAT);
            case 'FWE'
                xSPM.u  = spm_uc(pThr{1}, xSPM.df, xSPM.STAT, xSPM.R, 1, xSPM.S);
            case 'FDR'
                xSPM.u  = spm_uc_FDR(pThr{1}, xSPM.df, xSPM.STAT, 1, VspmT, 0);
            otherwise
                error('Variable corr must be either "unc", "FWE", or "FDR".');
        end
        
        [x, y, z] = ind2sub(size(Tvals), find(Tvals > xSPM.u));
        XYZ = [x'; y'; z'];
        Z = min(Inf, spm_get_data(SPM.xCon(iCon).Vspm, XYZ));

        V2R = 1/prod(SPM.xVol.FWHM(SPM.xVol.DIM > 1));
        [uc, xSPM.Pc, ue] = spm_uc_clusterFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, Z, XYZ, V2R, xSPM.u);
        [up, xSPM.Pp] = spm_uc_peakFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, Z, XYZ, xSPM.u);

        A = spm_clusters(XYZ);
        Q = [];
        for i = 1:max(A)
            j = find(A == i);
            if length(j) >= xSPM.k
                Q = [Q j];
            end
        end

        xSPM.Z     = Z(:,Q);
        xSPM.XYZ   = XYZ(:,Q);
        xSPM.XYZmm = xSPM.XYZ;

        for jVox = 1:length(xSPM.XYZ)
            xSPM.XYZmm(:, jVox) = (xSPM.XYZ(:, jVox).' * xSPM.M(1:3, 1:3) + xSPM.M(1:3, 4)').';
        end
        uu = spm_uc(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.S);

        xSPM.uc = [uu up ue uc];
        
        table = spm_list_edited(xSPM, 5, 8);

        if length(corr) == 2
            switch corr{2}
                case 'unc'
                    col = 4;
                case 'FWE'
                    col = 1;
                case 'FDR'
                    col = 2;
            end
        end

        % Cluster row index.
        clustLevelCol = table(:, 1);
        clustIdx = find(~cellfun(@isempty, clustLevelCol));

        % Pare down table to cluster-level significant clusters if
        % cluster-extent based thresholding is being applied. Else, keep as
        % is.
        if length(pThr) == 2
            for iClust = length(clustIdx):-1:1
                if iClust == length(clustIdx)
                    if table{clustIdx(iClust), col} > pThr{2}
                        table(clustIdx(iClust):end, :) = [];
                    end
                else
                    if table{clustIdx(iClust), col} > pThr{2}
                        table(clustIdx(iClust):clustIdx(iClust + 1) - 1, :) = [];
                    end
                end
            end
        end
        clustLevelCol = table(:, 1);
        clustIdx = find(~cellfun(@isempty, clustLevelCol));
                
        clear x y z XYZ
        
        clustHeader = VspmT;
        rawData = spm_read_vols(VspmT);
        [allClustVals, nClusters] = spm_bwlabel(double(rawData > xSPM.u), 18);

        [clustSize, clustNum] = sort(histc(allClustVals(:), 0:nClusters), 1, 'descend');
        clustSize = clustSize(2:end); clustNum = clustNum(2:end) - 1;
        clustNum = clustNum(clustSize >= k); clustSize = clustSize(clustSize >= k);

        for jClust = 1:length(clustNum)
            oneClustVals = zeros(size(allClustVals));
            oneClustVals(allClustVals == clustNum(jClust)) = 1;
            [x, y, z] = ind2sub(size(oneClustVals), find(oneClustVals == 1));
            clustXYZ = [x'; y'; z'];
            adjClustXYZ = adjust_XYZ(clustXYZ, clustHeader.mat, dHeader);
            clustNumber = 0;
            
            % Determine which cluster you're looking at and fill in that
            % position in the output table. Basically we're merging the
            % information about the clusters from the table with the
            % information found through more direct means, including effect
            % size.
            for kClust = 1:length(clustIdx)
                peakMM(1:3) = (table{clustIdx(kClust), 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);
                if sum(ismember(adjClustXYZ{1}.', peakMM, 'rows'))
                    clustNumber = kClust;
                    clustPeakMM = peakMM;
                end
            end
            if clustNumber ~= 0
                peakCoord = table{clustIdx(clustNumber), 10}.';
                clustHeader.fname = [outDir 'Cluster_' sprintf('%03d', clustNumber) '_' num2str(peakCoord(1)) '_' num2str(peakCoord(2)) '_' num2str(peakCoord(3)) '.nii'];
                spm_write_vol(clustHeader, oneClustVals);
                
                % Fill in output csv.
                outStruct{roaCol}.col{clustIdx(clustNumber), 1} = '';
                outStruct{baCol}.col{clustIdx(clustNumber), 1} = '';
                
                if table{clustIdx(clustNumber), 10}(1) < 0
                    outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'L';
                elseif table{clustIdx(clustNumber), 10}(1) > 0
                    outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'R';
                else
                    outStruct{lrCol}.col{clustIdx(clustNumber), 1} = 'I';
                end
                
                outStruct{kCol}.col{clustIdx(clustNumber), 1} = clustSize(jClust) * voxelScalar;
                outStruct{tCol}.col{clustIdx(clustNumber), 1} = num2str(table{clustIdx(clustNumber), 7});
                outStruct{dCol}.col{clustIdx(clustNumber), 1} = Dvals(clustPeakMM(1), clustPeakMM(2), clustPeakMM(3));
                outStruct{xCol}.col{clustIdx(clustNumber), 1} = peakCoord(1);
                outStruct{yCol}.col{clustIdx(clustNumber), 1} = peakCoord(2);
                outStruct{zCol}.col{clustIdx(clustNumber), 1} = peakCoord(3);
                
                if clustNumber == length(clustIdx)
                    [nRows, ~] = size(table);
                else
                    nRows = clustIdx(clustNumber + 1) - 1;
                end
                
                for mSubClust = clustIdx(clustNumber) + 1:nRows
                    outStruct{roaCol}.col{mSubClust, 1} = '';
                    outStruct{baCol}.col{mSubClust, 1} = '';
                    
                    if table{mSubClust, 10}(1) < 0
                        outStruct{lrCol}.col{mSubClust, 1} = 'L';
                    elseif table{mSubClust, 10}(1) > 0
                        outStruct{lrCol}.col{mSubClust, 1} = 'R';
                    else
                        outStruct{lrCol}.col{mSubClust, 1} = 'I';
                    end
                    
                    outStruct{kCol}.col{mSubClust, 1} = '';
                    outStruct{tCol}.col{mSubClust, 1} = num2str(table{mSubClust, 7});
                    subPeakMM = (table{mSubClust, 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);
                    outStruct{dCol}.col{mSubClust, 1} = Dvals(subPeakMM(1), subPeakMM(2), subPeakMM(3));
                    outStruct{xCol}.col{mSubClust, 1} = table{mSubClust, 10}(1);
                    outStruct{yCol}.col{mSubClust, 1} = table{mSubClust, 10}(2);
                    outStruct{zCol}.col{mSubClust, 1} = table{mSubClust, 10}(3);
                end
            end
        end

        fprintf('\t\t%d out of %d clusters are larger than %d voxels.\n', length(clustIdx), nClusters, k);
        write_csv(outStruct, [outDir 'ClusterReport_' sprintf('%03d', length(clustIdx)) '.csv']);
    else
        fprintf(['\tSkipping contrast ' num2str(iCon) ', ' SPM.xCon(iCon).name ', because it is an F con.\n']);
    end
end
cd(origDir);
fprintf('Done.\n\n');
end

%% Determine Voxel Size in Nifti File
function [sizeX, sizeY, sizeZ] = get_voxel_size(niiFile)
% FORMAT [sizeX, sizeY, sizeZ] = get_voxel_size(niiFile)
% Determines size of voxel in X, Y, and Z dimensions. This assumes that the
% matrix is not diagonal, but works fine if it is.
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

%% Adjust Coordinates of ROI
function funcXYZ = adjust_XYZ(XYZ, ROImat, V)
% FORMAT funcXYZ = adjust_XYZ(XYZ, ROImat, V)
% By Dennis Thompson.
%
%
% XYZ:              XYZ coordinates of ones in binary matrix.
% ROImat:           Transformation matrix from ROI file.
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
