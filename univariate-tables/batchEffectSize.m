% Example wrapper to call save_clusters_and_effect_size.m function.
spmFile = '/nfs/ep2/AX/analysis/MTU/func_an_SPM8/150406_SOBP/EP32vsHC32_paired_t-test/CueB-CueA/SPM.mat';

maskFile = '';

pThr = {0.005};
correction = {'unc'};
k = 10; % minimum cluster size

saveClustersAndEffectSize(spmFile, pThr, correction, k, maskFile);
