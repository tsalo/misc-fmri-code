% batch_newLSS
% General settings
tStart = tic;
subjects = {'epc141'};

spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
spmSubFolder = '/MTU/func_an_SPM8/';

outFolder = '/home/tsalo/newLSS/';
outSubFolder = '/';

% LSS settings
ignoreConditions = {'ProbeAX' 'ProbeAY' 'ProbeBX' 'ProbeBY' 'ExcludedCues' 'ExcludedProbes'};
settings.estimate = 1;         % 0- do not estimate SPM models, 1- estimate
settings.model = 2;            % 1- Rissman, 2- LSS
settings.deleteFiles = 0;      % 0- do not delete files, 1- delete files
settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings
rois = {'/nfs/ep2/masks/072313_EPC_BP_SZ/L_Fusiform_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/L_Inf_Parietal_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/gtmasks/I_Ant_Cingulum_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/gtmasks/L_Mid_Occipital_r5.nii'
        '/nfs/ep2/masks/072313_EPC_BP_SZ/gtmasks/R_Mid_Occipital_r5.nii'};
settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi

for iSubj = 1:length(subjects)
    spmDir = [spmFolder subjects{iSubj} spmSubFolder];
    outDir = [outFolder subjects{iSubj} outSubFolder];
    
    images = generate_spm_singletrial_newLSS(subjects{iSubj}, spmDir, outDir, ignoreConditions, settings);
    newLSS_correlation(images, [spmDir '/SPM.mat'], rois, settings);
end
hourToc(tStart);