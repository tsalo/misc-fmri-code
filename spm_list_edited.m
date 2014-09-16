function outTable = spm_list_edited(xSPM, Num, Dis)
% FORMAT outTable = spm_list_edited(xSPM, Num, Dis)
% Summary list of local maxima for entire volume of interest
%
% Inputs:
% xSPM   - structure containing SPM, distribution & filtering details
%        - required fields are:
%    .Z     - minimum of n Statistics {filtered on u and k}
%    .n     - number of conjoint tests CAN DO IGNORE CONJUNCTIONS AND SET TO 1
%    .STAT  - distribution {Z, T, X or F} CAN DO
%    .df    - degrees of freedom [df{interest}, df{residual}] CAN DO
%    .u     - height threshold {Z} CAN DO
%    .k     - extent threshold {voxels} CAN DO (preset)
%    .XYZ   - location of voxels {voxel coords} CAN DO (same as ROIs)
%    .XYZmm - location of voxels {mm coords} CAN DO (convert XYZ to mm using matrix)
%    .S     - search Volume {voxels}
%    .R     - search Volume {resels}
%    .FWHM  - smoothness {voxels}
%    .M     - voxels - > mm matrix CAN DO
%    .VOX   - voxel dimensions {mm} CAN DO
%    .DIM   - image dimensions {voxels} CAN DO
%    .units - space units CAN DO SET TO {'mm' 'mm' 'mm'}
%    .VRpv  - filehandle - Resels per voxel CAN DO (load RPV.img header)
%    .Ps    - uncorrected P values in searched volume (for voxel FDR)
%    .Pp    - uncorrected P values of peaks (for peak FDR) 
%    .Pc    - uncorrected P values of cluster extents (for cluster FDR)
%    .uc    - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc MAY NOT NEED
%    .thresDesc - description of height threshold (string) DO NOT NEED
%
% Num    - number of maxima per cluster [3]
% Dis    - distance among clusters {mm} [8]
%
% Output:
% outTable  - Nx9 cell array (same as whole brain results from the GUI
%             without the set level stats).
%__________________________________________________________________________
%
% spm_list_edited characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected subsets
% or clusters).  The corrected significance of the results are based on
% set, cluster and voxel-level inferences using distributional
% approximations from the Theory of Gaussian Fields.  These distributions
% assume that the SPM is a reasonable lattice approximation of a
% continuous random field with known component field smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%               - Qc   = lowest FDR bound for which this cluster would be
%                        declared positive
%
% Peak-level    - T/F  = Statistic upon which the SPM is based
%               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qp   = lowest FDR bound for which this peak would be
%                        declared positive
%               - Pu   = prob(Ze or higher at that voxel)
%
% Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the MATLAB
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the MATLAB workspace.
%
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Andrew Holmes, Guillaume Flandin
% $Id: spm_list.m 4617 2012-01-11 15:46:16Z will $
% Edited and reduced by Taylor Salo 140912 for compatibility with
% save_clusters_and_effect_size.

%-Get number of maxima per cluster to be reported
%----------------------------------------------------------------------
if ~exist('Num', 'var')
    Num = spm_get_defaults('stats.results.volume.nbmax'); % 3
end

%-Get minimum distance among clusters (mm) to be reported
%----------------------------------------------------------------------
if ~exist('Dis', 'var')
    Dis = spm_get_defaults('stats.results.volume.nbmax'); % 8
end

%-Get header string
%----------------------------------------------------------------------
if xSPM.STAT ~= 'P'
    Title = 'p-values adjusted for search volume';
else
    Title = 'Posterior Probabilities';
end

%-Extract data from xSPM
%----------------------------------------------------------------------
S         = xSPM.S;
VOX       = xSPM.VOX;
DIM       = xSPM.DIM;
M         = xSPM.M;
XYZ       = xSPM.XYZ;
Z         = xSPM.Z;
VRpv      = xSPM.VRpv;
n         = xSPM.n;
STAT      = xSPM.STAT;
df        = xSPM.df;
u         = xSPM.u;
k         = xSPM.k;
try, uc   = xSPM.uc; end
try, QPs  = xSPM.Ps; end
try, QPp  = xSPM.Pp; end
try, QPc  = xSPM.Pc; end

if STAT~='P'
    R     = full(xSPM.R);
    FWHM  = full(xSPM.FWHM);
end
try
    units = xSPM.units;
catch
    units = {'mm' 'mm' 'mm'};
end
units{1}  = [units{1} ' '];
units{2}  = [units{2} ' '];

DIM       = DIM > 1;              % non-empty dimensions
D         = sum(DIM);             % highest dimension
VOX       = VOX(DIM);             % scaling

if STAT ~= 'P'
    FWHM  = FWHM(DIM);            % Full width at half max
    FWmm  = FWHM.*VOX;            % FWHM {units}
    V2R   = 1/prod(FWHM);         % voxels to resels
    k     = k*V2R;                % extent threshold in resels
    R     = R(1:(D + 1));         % eliminate null resel counts
    try, QPs = sort(QPs(:)); end  % Needed for voxel   FDR
    try, QPp = sort(QPp(:)); end  % Needed for peak    FDR
    try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
end

% Choose between voxel-wise and topological FDR
%----------------------------------------------------------------------
topoFDR = spm_get_defaults('stats.topoFDR');

%-Tolerance for p-value underflow, when computing equivalent Z's
%----------------------------------------------------------------------
tol = eps*10;

%-Table Headers
%----------------------------------------------------------------------
TabDat.tit = Title;

TabDat.hdr = {...
    'set',      'p',            '\itp';...
    'set',      'c',            '\itc';...
    'cluster',  'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
    'cluster',  'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
    'cluster',  'equivk',       '\itk\rm_E';...
    'cluster',  'p(unc)',       '\itp\rm_{uncorr}';...
    'peak',     'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
    'peak',     'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
    'peak',      STAT,          sprintf('\\it%c',STAT);...
    'peak',     'equivZ',       '(\itZ\rm_\equiv)';...
    'peak',     'p(unc)',       '\itp\rm_{uncorr}';...
    '',         'x,y,z {mm}',   [units{:}]}';...

%-Coordinate Precisions
%----------------------------------------------------------------------
if isempty(xSPM.XYZmm) % empty results
    xyzfmt = '%3.0f %3.0f %3.0f';
    voxfmt = repmat('%0.1f ',1,numel(FWmm));
elseif ~any(strcmp(units{3},{'mm',''})) % 2D data
    xyzfmt = '%3.0f %3.0f %3.0f';
    voxfmt = repmat('%0.1f ',1,numel(FWmm));
else % 3D data, work out best precision based on voxel sizes and FOV
    xyzsgn = min(xSPM.XYZmm(DIM,:),[],2) < 0;
    xyzexp = max(floor(log10(max(abs(xSPM.XYZmm(DIM,:)),[],2)))+(max(abs(xSPM.XYZmm(DIM,:)),[],2) >= 1),0);
    voxexp = floor(log10(abs(VOX')))+(abs(VOX') >= 1);
    xyzdec = max(-voxexp,0);
    voxdec = max(-voxexp,1);
    xyzwdt = xyzsgn+xyzexp+(xyzdec>0)+xyzdec;
    voxwdt = max(voxexp,0)+(voxdec>0)+voxdec;
    tmpfmt = cell(size(xyzwdt));
    for i = 1:numel(xyzwdt)
        tmpfmt{i} = sprintf('%%%d.%df ', xyzwdt(i), xyzdec(i));
    end
    xyzfmt = [tmpfmt{:}];
    tmpfmt = cell(size(voxwdt));
    for i = 1:numel(voxwdt)
        tmpfmt{i} = sprintf('%%%d.%df ', voxwdt(i), voxdec(i));
    end
    voxfmt = [tmpfmt{:}];
end
TabDat.fmt = {  '%-0.3f','%g',...                          %-Set
    '%0.3f', '%0.3f','%0.0f', '%0.3f',...                  %-Cluster
    '%0.3f', '%0.3f', '%6.2f', '%5.2f', '%0.3f',...        %-Peak
    xyzfmt};                                               %-XYZ

%-Table filtering note
%----------------------------------------------------------------------
if isinf(Num)
    TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
else
    TabDat.str = sprintf(['table shows %d local maxima ',...
        'more than %.1fmm apart'],Num,Dis);
end 

%-Footnote with SPM parameters
%----------------------------------------------------------------------
if STAT ~= 'P'
    Pz              = spm_P(1,0,u,df,STAT,1,n,S);
    Pu              = spm_P(1,0,u,df,STAT,R,n,S);
    [P Pn Ec Ek]    = spm_P(1,k,u,df,STAT,R,n,S);

    TabDat.ftr      = cell(9,2);
    TabDat.ftr{1,1} = ...
        ['Height threshold: ' STAT ' = %0.2f, p = %0.3f (%0.3f)'];
    TabDat.ftr{1,2} = [u,Pz,Pu];
    TabDat.ftr{2,1} = ...
        'Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)';
    TabDat.ftr{2,2} = [k/V2R,Pn,P];
    TabDat.ftr{3,1} = ...
        'Expected voxels per cluster, <k> = %0.3f';
    TabDat.ftr{3,2} = Ek/V2R;
    TabDat.ftr{4,1} = ...
        'Expected number of clusters, <c> = %0.2f';
    TabDat.ftr{4,2} = Ec*Pn;
    if any(isnan(uc))
        TabDat.ftr{5,1} = 'FWEp: %0.3f, FDRp: %0.3f';
        TabDat.ftr{5,2} = uc(1:2);
    else
        TabDat.ftr{5,1} = ...
            'FWEp: %0.3f, FDRp: %0.3f, FWEc: %0.0f, FDRc: %0.0f';
        TabDat.ftr{5,2} = uc;
    end
    TabDat.ftr{6,1} = 'Degrees of freedom = [%0.1f, %0.1f]';
    TabDat.ftr{6,2} = df;
    TabDat.ftr{7,1} = ...
        ['FWHM = ' voxfmt units{:} '; ' voxfmt '{voxels}'];
    TabDat.ftr{7,2} = [FWmm FWHM];
    TabDat.ftr{8,1} = ...
        'Volume: %0.0f = %0.0f voxels = %0.1f resels';
    TabDat.ftr{8,2} = [S*prod(VOX),S,R(end)];
    TabDat.ftr{9,1} = ...
        ['Voxel size: ' voxfmt units{:} '; (resel = %0.2f voxels)'];
    TabDat.ftr{9,2} = [VOX,prod(FWHM)];
 else
    TabDat.ftr = {};
end 

%-Characterize excursion set in terms of maxima
% (sorted on Z values and grouped by regions)
%----------------------------------------------------------------------
if isempty(Z)
    TabDat.dat = cell(0,12);
    outTable = TabDat.dat(:, 3:end);
    return
end

%-Workaround in spm_max for conjunctions with negative thresholds
%----------------------------------------------------------------------
minz           = abs(min(min(Z)));
Z              = 1 + minz + Z;
[N Z XYZ A L]  = spm_max(Z,XYZ);
Z              = Z - minz - 1;

%-Convert cluster sizes from voxels (N) to resels (K)
%----------------------------------------------------------------------
c              = max(A);                           %-Number of clusters
NONSTAT        = spm_get_defaults('stats.rft.nonstat');

if STAT ~= 'P'
    if NONSTAT
        K      = zeros(c,1);
        for i  = 1:c

            %-Get LKC for voxels in i-th region
            %----------------------------------------------------------
            LKC = spm_get_data(VRpv,L{i});

            %-Compute average of valid LKC measures for i-th region
            %----------------------------------------------------------
            valid = ~isnan(LKC);
            if any(valid)
                LKC = sum(LKC(valid)) / sum(valid);
            else
                LKC = V2R; % fall back to whole-brain resel density
            end

            %-Intrinsic volume (with surface correction)
            %----------------------------------------------------------
            IV   = spm_resels([1 1 1],L{i},'V');
            IV   = IV*[1/2 2/3 2/3 1]';
            K(i) = IV*LKC;

        end
        K = K(A);
    else
        K = N*V2R;
    end
end

%-Convert maxima locations from voxels to mm
%----------------------------------------------------------------------
XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];

%-Set-level p values {c} - do not display if reporting a single cluster
%----------------------------------------------------------------------
if STAT ~= 'P'
    Pc     = spm_P(c,k,u,df,STAT,R,n,S);            %-Set-level p-value
else
    Pc     = [];
end

TabDat.dat = {Pc,c};
TabLin     = 1;

%-Cluster and local maxima p-values & statistics
%----------------------------------------------------------------------
while numel(find(isfinite(Z)))

    %-Find largest remaining local maximum
    %------------------------------------------------------------------
    [U,i]  = max(Z);            %-largest maxima
    j      = find(A == A(i));   %-maxima in cluster


    %-Compute cluster {k} and peak-level {u} p values for this cluster
    %------------------------------------------------------------------
    if STAT ~= 'P'

        % p-values (FWE)
        %--------------------------------------------------------------
        Pz      = spm_P(1,0,   U,df,STAT,1,n,S);  % uncorrected p value
        Pu      = spm_P(1,0,   U,df,STAT,R,n,S);  % FWE-corrected {based on Z}
        [Pk Pn] = spm_P(1,K(i),u,df,STAT,R,n,S);  % [un]corrected {based on K}

        % q-values (FDR)
        %--------------------------------------------------------------
        if topoFDR
            Qc  = spm_P_clusterFDR(K(i),df,STAT,R,n,u,QPc); % based on K
            Qp  = spm_P_peakFDR(U,df,STAT,R,n,u,QPp);       % based on Z
            Qu  = [];
        else
            Qu  = spm_P_FDR(U,df,STAT,n,QPs);     % voxel FDR-corrected
            Qc  = [];
            Qp  = [];
        end

        % Equivalent Z-variate
        %--------------------------------------------------------------
        if Pz < tol
            Ze  = Inf;
        else
            Ze  = spm_invNcdf(1 - Pz);
        end
    else
        Pz      = [];
        Pu      = [];
        Qu      = [];
        Pk      = [];
        Pn      = [];
        Qc      = [];
        Qp      = [];
        ws      = warning('off','SPM:outOfRangeNormal');
        Ze      = spm_invNcdf(U);
        warning(ws);
    end

    if topoFDR
    [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu,Qp,U,Ze,Pz,XYZmm(:,i));
    else
    [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu,Qu,U,Ze,Pz,XYZmm(:,i));
    end
    TabLin = TabLin + 1;

    %-Print Num secondary maxima (> Dis mm apart)
    %------------------------------------------------------------------
    [~,q] = sort(-Z(j));                              % sort on Z value
    D     = i;
    for i = 1:length(q)
        d = j(q(i));
        if min(sqrt(sum((XYZmm(:,D)-repmat(XYZmm(:,d),1,size(D,2))).^2)))>Dis
            if length(D) < Num
                % voxel-level p values {Z}
                %------------------------------------------------------
                if STAT ~= 'P'
                    Pz     = spm_P(1,0,Z(d),df,STAT,1,n,S);
                    Pu     = spm_P(1,0,Z(d),df,STAT,R,n,S);
                    if topoFDR
                        Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
                        Qu = [];
                    else
                        Qu = spm_P_FDR(Z(d),df,STAT,n,QPs);
                        Qp = [];
                    end
                    if Pz < tol
                        Ze = Inf;
                    else
                        Ze = spm_invNcdf(1 - Pz); 
                    end
                else
                    Pz     = [];
                    Pu     = [];
                    Qu     = [];
                    Qp     = [];
                    ws     = warning('off','SPM:outOfRangeNormal');
                    Ze     = spm_invNcdf(Z(d));
                    warning(ws);
                end
                D     = [D d];
                if topoFDR
                [TabDat.dat{TabLin,7:12}] = ...
                    deal(Pu,Qp,Z(d),Ze,Pz,XYZmm(:,d));
                else
                [TabDat.dat{TabLin,7:12}] = ...
                    deal(Pu,Qu,Z(d),Ze,Pz,XYZmm(:,d));
                end
                TabLin = TabLin+1;
            end
        end
    end
    Z(j) = NaN;     % Set local maxima to NaN
end
outTable = TabDat.dat(:, 3:end);
end
