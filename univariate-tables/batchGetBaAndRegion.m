% batchGetBaAndRegion
clusterImage = '/nfs/ep2/AX/analysis/MTU/func_an_SPM8/021113_HC30_SZ31/HCvsSZ/CueB-CueA_unequalvariance/0.001_unc_clusters/SZ-HC/Cluster_001.nii';

[outString, baString] = getBaAndRegion(clusterImage);
