# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:40:54 2014
Summarizes betas and calls extract_betas.
@author: taylorsalo
"""

import edit_csv as ec
import extract_betas as eb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

subject_list = "/nfs/ep2/AX/project_folders/072313_EPC_BP_SZ/subject_list.csv"
file_sup_path = "/nfs/ep2/AX/first_levels/00_MONTH/"
file_sub_path = "/MTU/func_an_SPM8/"
mask_file = "/nfs/ep2/masks/PickAtlas/B46_L_d0.nii"
group_list = ["HC", "SZ", "BP"]

subjects = ec.read_csv(subject_list, group_list)

wanted_contrasts = eb.get_spm_contrast_list(file_sup_path + subjects[0][1] +
                                            file_sub_path + "SPM.mat")

headers = ["ID", "Group"] + wanted_contrasts
out_struct = pd.DataFrame(columns=headers)
for group in subjects:
    for j_subj in range(1, len(group)):
        beta_values = []
        wanted_file_list = eb.get_spm_contrast_files(file_sup_path +
                                                     group[j_subj] +
                                                     file_sub_path +
                                                     "SPM.mat",
                                                     wanted_contrasts)
        for k_con, contrast in enumerate(wanted_contrasts):
            beta_file = (file_sup_path + group[j_subj] + file_sub_path +
                         wanted_file_list[k_con])
            print(beta_file)
            beta_value = eb.main(beta_file, mask_file)
            beta_values.append(beta_value)
        out_struct.loc[len(out_struct)] = ([group[j_subj], group[0]] +
                                           beta_values)

grouped = out_struct.groupby('Group')
for con in wanted_contrasts:
    out_struct[con] = out_struct[con].astype(float)

mean_ = grouped.aggregate(lambda x: np.mean(x)).transpose()
mean_ = mean_.reindex_axis(group_list, axis=1)
mean_ = mean_.reindex_axis(wanted_contrasts, axis=0)
ste_ = grouped.aggregate(lambda x: np.std(x) / np.sqrt(x.count())).transpose()
ste_ = ste_.reindex_axis(group_list, axis=1)
ste_ = ste_.reindex_axis(wanted_contrasts, axis=0)

max_value = (mean_ + ste_).max().max()
min_value = (mean_ - ste_).min().min()
range_ = max_value - min_value

# Use max_value, min_value, and range_ to determine scale for y-axis.


ax = mean_.plot(kind='bar', yerr=ste_, rot=45)
plt.axhline(0, color='k')
plt.tight_layout()
