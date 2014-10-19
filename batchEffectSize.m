% Example wrapper to call save_clusters_and_effect_size.m function.
spmFile = '/nfs/cntracs/analysis/MTC/func_an_SPM8_bestproc/070814_paired_t-test/HC/CueB-CueA/SPM.mat';

pThr = {0.01 0.05};
corr = {'unc' 'FWE'};
k = 0;

saveClustersAndEffectSize(spmFile, pThr, corr, k);
