# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:40:54 2014
Summarizes betas and calls extract_betas.
@author: taylorsalo
"""

import csv
import extract_betas as eb
import pandas as pd


def transpose_list_of_lists(list_):
    transposed_list_ = [[row[col] for row in list_ if row[col]] for col in
                        range(len(list_[0]))]
    transposed_list = [col for col in transposed_list_ if col]
    return transposed_list

subject_list = "/nfs/ep2/AX/project_folders/072313_EPC_BP_SZ/subject_list.csv"
file_sup_path = "/nfs/ep2/AX/first_levels/00_MONTH/"
file_sub_path = "/MTU/func_an_SPM8/"
file_name = "con_0001.img"
mask_file = "/nfs/ep2/masks/PickAtlas/B46_L_d0.nii"
groups = ["HC", "SZ", "BP"]

subjects = []

with open(subject_list, 'r') as fo:
    reader = csv.reader(fo)
    for row in reader:
        subjects.append(row)

subjects = transpose_list_of_lists(subjects)

wanted_contrasts = eb.determine_spm_contrasts(file_sup_path + subjects[0][1] +
                                              file_sub_path + "SPM.mat")

headers = ["ID", "Group"] + wanted_contrasts
out_struct = pd.DataFrame(columns=headers)
for i_group in subjects:
    if i_group[0] in groups:
        for j in range(1, len(i_group)):
            
            beta_file = file_sup_path + i_group[j] + file_sub_path + file_name
            beta_value = eb.main(beta_file, mask_file)
            betas = [beta_value]
            out_struct.loc[len(out_struct)] = [i_group[j], i_group[0]] + betas
