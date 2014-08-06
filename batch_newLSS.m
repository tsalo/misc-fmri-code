% batch_newLSS
% General settings
tStart = tic;
subjects = {'epc141'};

spmDir = '/nfs/ep2/AX/first_levels/00_MONTH/';
spmSubDir = '/MTU/func_an_SPM8/';

outDir = '/home/tsalo/newLSS/';
outSubDir = '/';

% LSS settings
excluded_conditions = {'ProbeAX' 'ProbeAY' 'ProbeBX' 'ProbeBY' 'ExcludedCues' 'ExcludedProbes'};
settings.estimate = 1;         % 0- do not estimate SPM models, 1- estimate
settings.modeltype = 2;        % 1- Rissman, 2- LSS
settings.discard_mm_files = 0; % 0- do not delete files, 1- delete files
settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings
rois = {'/nfs/ep2/masks/072313_EPC_BP_SZ/L_Fusiform_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/L_Inf_Parietal_r5.nii'};
settings.fconnType = 'seed2voxel'; % seed2voxel or roi2roi

for iSubj = 1:length(subjects)
    spmFolder = [spmDir subjects{iSubj} spmSubDir];
    outFolder = [outDir subjects{iSubj} outSubDir];
    
    images = generate_spm_singletrial_newLSS(subjects{iSubj}, spmFolder, outFolder, excluded_conditions, settings);
    newLSS_correlation(images, [spmFolder '/SPM.mat'], rois, settings);
end

hourToc(tStart);