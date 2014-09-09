Misc_fMRI_code
==============

Just small functions for preprocessing, analysis, etc. of fMRI data.

* Batch_newLSS/generate_spm_singletrial_newLSS/newLSS_correlation: Scripts to perform the Least Squares- Separate procedure from Mumford et al. 2012 and improved in Turner et al. 2012, as well as to perform functional connectivity analysis on the results. Generate_spm_singletrial, written by Maureen Ritchey, matched the earlier version of LSS, and I altered it to match the later version. The functional connectivity portion (newLSS_correlation) is largely an amalgamation of several scripts written by Dennis Thompson, with some slight alterations to make it compatible with LSS. Batch_newLSS is just a wrapper that sets paths and settings and calls both generate_spm_singletrial_newLSS and newLSS_correlation.
  * To do:
    1. Figure out graph theory analyses to apply to roi2roi results.
* Clustsim calls AFNI functions to perform Monte Carlo simulations on a second-level SPM analysis. The paths to the AFNI functions will need to be altered by anyone using the function.
  * To do:
    1. Add in code to read the outputted text file so that the results can be outputted as a variable. 
* Save_clusters_and_effect_size takes a second-level analysis, creates a Cohen’s d image for each t contrast, creates masks of each significant cluster in each t contrast, and extracts the mean Cohen’s d for each cluster, summarizing it all in csv file. The Cohen’s d image can also be used to determine the Cohen’s d of the peak voxel in each cluster.
  * To do:
    1. Add in cluster level stats for cluster-extent thresholding.
    2. Reorganize outputted csv tables so that they match our lab's results tables.
    3. Include subpeaks?
