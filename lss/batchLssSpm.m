% batch_newLSS
% General settings
tStart = tic;
subjects = {'epc138'};

spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
spmSubFolder = '/MTU/func_an_SPM8/';

outFolder = '/home/tsalo/lssForUnivariate/';
outSubFolder = '/';

% LSS settings
includeConditions = {'CueA' 'CueB'};
settings.model = 2;            % 1- Rissman, 2- LSS
settings.useTempFS = 1;        % 0- do not use temporary files, 1- use temporary files
settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings
rois = load('/nfs/cntracs/lssForKim/gt_rois_final.mat');
settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi

for iSubj = 1:length(subjects)
    spmDir = [spmFolder subjects{iSubj} spmSubFolder];
    outDir = [outFolder subjects{iSubj} outSubFolder];
    
    images = lssGenerateBetasSpm(subjects{iSubj}, spmDir, outDir, includeConditions, settings);
    lssCorrelation(images, rois.rois, settings);
end
hourToc(tStart);