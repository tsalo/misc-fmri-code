% Example wrapper to call save_clusters_and_effect_size.m function.
spmFile = '/nfs/ep/EP_Processing/AX/analysis/err_MTU/func_an_SPM8/070214_AXonly_lapse_SZ/SZ32HC23_AXlapse/SZ/AX_CuCo_PrIn_Cu-AX_CuCo_PrCo_Cu/SPM.mat';

pThr = {0.01 0.05};
corr = {'unc' 'FDR'};
k = 0;

saveClustersAndEffectSize(spmFile, pThr, corr, k);