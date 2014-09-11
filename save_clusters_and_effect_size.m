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
% pThr:    p threshold. Cell array (1x1 or 1x2) of doubles.
% corr:    Correction applied to p threshold (FWE, FDR, or unc). Cell array
%          (same length as pThr) of strings.
% k:       Minimum acceptable cluster size. Double.

if exist(spmFile, 'file')
    [path, ~] = fileparts(spmFile);
    load(spmFile);
else
    error('spmFile does not exist. Quitting.');
end

% if iscell(pThr)
%     fprintf('Warning, your set p threshold is a string. Trying to convert to double.\n')
%     pThr = str2double(pThr);
%     if isnan(pThr)
%         fprintf('P threshold could not be converted to double, setting to default 0.001.\n');
%         pThr = 0.001;
%     end
% elseif ~isa(pThr, 'double')
%     fprintf('Warning, your set p threshold is not a double or a string. Setting to default 0.001.\n');
%     pThr = 0.001;
% end
% 
% if ischar(k)
%     fprintf('Warning, your set k is a string. Trying to convert to double.\n')
%     k = str2double(k);
%     if isnan(k)
%         fprintf('P threshold could not be converted to double, setting to default 5.\n');
%         k = 5;
%     end
% elseif ~isa(k, 'double')
%     fprintf('Warning, your set k is not a double or a string. Setting to default 5.\n');
%     k = 5;
% end

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
        
        spmT = [path '/' SPM.xCon(iCon).Vspm.fname];
        conName = ['Contrast_' sprintf('%03d', iCon) '-' strrep(SPM.xCon(iCon).name, ' ', '_')];
        fprintf(['\tEvaluating contrast ' num2str(iCon) ', ' conName '\n']);
        
        if length(pThr) == 1
            outDir = [path '/v' num2str(pThr{1}) '_' corr{1} '_clusters/' conName '/'];
        else
            outDir = [path '/v' num2str(pThr{1}) '_' corr{1} '_c' num2str(pThr{2}) '_' corr{2} '_clusters/' conName '/'];
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
        xSPM.DIM = SPM.xVol.DIM;
        
        % Set up output csv.
        clear outStruct
        outStruct{1}.header{1} = 'Region of Activation'; outStruct{2}.header{1} = 'BA'; outStruct{3}.header{1} = 'L/R';
        outStruct{4}.header{1} = 'k (mm3)'; outStruct{5}.header{1} = 'T'; outStruct{6}.header{1} = 'D';
        outStruct{7}.header{1} = 'x'; outStruct{8}.header{1} = 'y'; outStruct{9}.header{1} = 'z';
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
        
        TabDat2 = spm_list_edited('List', xSPM);
        table = TabDat2.dat(:, 3:end);

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
        clustLevel = table(:, 1);
        clustIdx = find(~cellfun(@isempty, clustLevel));

        % Pare down table
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
        clustLevel = table(:, 1);
        clustIdx = find(~cellfun(@isempty, clustLevel));
                
        clear x y z XYZ
        
        clustHeader = VspmT;
        rawData = spm_read_vols(VspmT);
        [allClustMat, nClusters] = spm_bwlabel(double(rawData > xSPM.u), 18);

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
            clustNumber = 0;
            for kClust = 1:length(clustIdx)
                peakMM(1:3) = (table{clustIdx(kClust), 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);
                if sum(ismember(betaXYZ{1}.', peakMM, 'rows'))
                    clustNumber = kClust;
                    clustPeakMM = peakMM;
                end
            end
            if clustNumber ~= 0
                peakCoord = table{clustIdx(clustNumber), 10}.';
                clustHeader.fname = [outDir 'Cluster_' sprintf('%03d', clustNumber) '_' num2str(peakCoord(1)) '_' num2str(peakCoord(2)) '_' num2str(peakCoord(3)) '.nii'];
                spm_write_vol(clustHeader, oneClustMat);
                
                % Fill in output csv.
                outStruct{1}.col{clustIdx(clustNumber), 1} = '';
                outStruct{2}.col{clustIdx(clustNumber), 1} = '';
                
                
                if table{clustIdx(clustNumber), 10}(1) < 0
                    outStruct{3}.col{clustIdx(clustNumber), 1} = 'L';
                elseif table{clustIdx(clustNumber), 10}(1) > 0
                    outStruct{3}.col{clustIdx(clustNumber), 1} = 'R';
                else
                    outStruct{3}.col{clustIdx(clustNumber), 1} = 'I';
                end
                outStruct{4}.col{clustIdx(clustNumber), 1} = clustSize(jClust) * voxelScalar;
                outStruct{5}.col{clustIdx(clustNumber), 1} = num2str(table{clustIdx(clustNumber), 7});
                outStruct{6}.col{clustIdx(clustNumber), 1} = Dvals(clustPeakMM(1), clustPeakMM(2), clustPeakMM(3));
                outStruct{7}.col{clustIdx(clustNumber), 1} = peakCoord(1);
                outStruct{8}.col{clustIdx(clustNumber), 1} = peakCoord(2);
                outStruct{9}.col{clustIdx(clustNumber), 1} = peakCoord(3);
                
                if clustNumber == length(clustIdx)
                    [nRows, ~] = size(table);
                else
                    nRows = clustIdx(clustNumber + 1) - 1;
                end
                for mSubClust = clustIdx(clustNumber) + 1:nRows
                    subPeakMM = (table{mSubClust, 10}.' - VspmT.mat(1:3, 4)') / VspmT.mat(1:3, 1:3);

                    outStruct{1}.col{mSubClust, 1} = '';
                    outStruct{2}.col{mSubClust, 1} = '';
                    if table{mSubClust, 10}(1) < 0
                        outStruct{3}.col{mSubClust, 1} = 'L';
                    elseif table{mSubClust, 10}(1) > 0
                        outStruct{3}.col{mSubClust, 1} = 'R';
                    else
                        outStruct{3}.col{mSubClust, 1} = 'I';
                    end
                    outStruct{4}.col{mSubClust, 1} = '';
                    outStruct{5}.col{mSubClust, 1} = num2str(table{mSubClust, 7});
                    outStruct{6}.col{mSubClust, 1} = Dvals(subPeakMM(1), subPeakMM(2), subPeakMM(3));
                    outStruct{7}.col{mSubClust, 1} = table{mSubClust, 10}(1);
                    outStruct{8}.col{mSubClust, 1} = table{mSubClust, 10}(2);
                    outStruct{9}.col{mSubClust, 1} = table{mSubClust, 10}(3);
                end
            end
        end

        fprintf('\t\t%d out of %d clusters are larger than %d voxels.\n', length(clustIdx), nClusters, k);
        write_csv(outStruct, [outDir 'ClusterReport_' sprintf('%03d', length(clustIdx)) '.csv']);
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
