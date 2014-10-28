% batch_newLSS
% General settings
tStart = tic;
subjects = {'epp343'};

spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
spmSubFolder = '/MTU/func_an_SPM8/';

outFolder = '/home/tsalo/lssForUnivariate/';
outSubFolder = '/';

% LSS settings
ignoreConditions = {'ProbeAX' 'ProbeAY' 'ProbeBX' 'ProbeBY' 'ExcludedCues' 'ExcludedProbes'};
settings.model = 2;            % 1- Rissman, 2- LSS
settings.deleteFiles = 0;      % 0- do not delete files, 1- delete files
settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings
rois = load('/nfs/cntracs/lssForKim/gt_rois2.mat');
settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi

for iSubj = 1:length(subjects)
    spmDir = [spmFolder subjects{iSubj} spmSubFolder];
    outDir = [outFolder subjects{iSubj} outSubFolder];
    
    images = lssGenerateBetasSpm(subjects{iSubj}, spmDir, outDir, ignoreConditions, settings);
%     lssCorrelation(images, rois.rois, settings);
end
hourToc(tStart);