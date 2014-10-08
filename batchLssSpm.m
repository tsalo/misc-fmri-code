% batch_newLSS
% General settings
tStart = tic;
subjects = {'epc141'};

spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
spmSubFolder = '/MTU/func_an_SPM8/';

outFolder = '/home/tsalo/newRissman/';
outSubFolder = '/';

% LSS settings
ignoreConditions = {'ProbeAX' 'ProbeAY' 'ProbeBX' 'ProbeBY' 'ExcludedCues' 'ExcludedProbes'};
settings.model = 1;            % 1- Rissman, 2- LSS
settings.deleteFiles = 0;      % 0- do not delete files, 1- delete files
settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings
rois = {'/nfs/ep2/masks/072313_EPC_BP_SZ/L_Fusiform_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/L_Inf_Parietal_r5.nii'};
settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi

for iSubj = 1:length(subjects)
    spmDir = [spmFolder subjects{iSubj} spmSubFolder];
    outDir = [outFolder subjects{iSubj} outSubFolder];
    
    images = lssGenerateBetasSpm(subjects{iSubj}, spmDir, outDir, ignoreConditions, settings);
    lssCorrelation(images, rois, settings);
end
hourToc(tStart);