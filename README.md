misc-fmri-code
==============

Just small functions for preprocessing, analysis, etc. of fMRI data.

* [batchLssSpm](batchLssSpm.m)/[lssGenerateBetasSpm](lssGenerateBetasSpm.m)/[lssCorrelation](lssCorrelation.m): Scripts to perform the Least Squares- Separate procedure from Mumford et al. 2012 and improved in Turner et al. 2012, as well as to perform functional connectivity analysis on the results. The base, [generate_spm_singletrial](https://github.com/ritcheym/fmri_misc/blob/master/generate_spm_singletrial.m), written by Maureen Ritchey, matched the earlier version of LSS, and I altered it to match the later version. The functional connectivity portion (lssCorrelation) is largely an amalgamation of several scripts written by Dennis Thompson, with some slight alterations to make it compatible with LSS. BatchLssSpm is just a wrapper that sets paths and settings and calls both lssGenerateBetasSpm and lssCorrelation.
  * To do:
    1. Figure out graph theory analyses to apply to roi2roi results.
* [clustsim](clustsim.m) calls AFNI functions to perform Monte Carlo simulations on a second-level SPM analysis. The paths to the AFNI functions will need to be altered by anyone using the function.
  * To do:
    1. Add in code to read the outputted text file so that the results can be outputted as a variable. 
* [batchEffectSize](univariate-tables/batchEffectSize.m)/[saveClustersAndEffectSize](univariate-tables/saveClustersAndEffectSize.m) take a second-level analysis, create a Cohen’s d image for each t contrast, create masks of each significant cluster in each t contrast, and extract the mean Cohen’s d for each cluster, summarizing it all in a csv file. The Cohen’s d image can also be used to determine the Cohen’s d of the peak voxel in each cluster. The script batchEffectSize is simply a wrapper to call the function saveClustersAndEffectSize. 
