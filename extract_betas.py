# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:39:30 2014
Extracts the mean value of the masked voxels from a nifti image.

@author: tsalo
"""

import sys
import os.path
import nibabel as nb
import numpy as np
import scipy.io


def main(nifti_image, mask_image):
    if os.path.isfile(nifti_image):
        nii_img = nb.load(nifti_image)
        nii_vals = np.array(nii_img.get_data())
        nii_header = nii_img.get_header()
        nii_shape = nii_header.get_data_shape()
        nii_affine = nii_header.get_qform()
    else:
        print("Variable nifti_image points to non-existent file.")
        return ""

    if os.path.isfile(mask_image):
        mask_img = nb.load(mask_image)
        mask_header = mask_img.get_header()
        mask_shape = mask_header.get_data_shape()
        mask_affine = mask_header.get_qform()

        if mask_shape == nii_shape:
            mask_index = np.where(mask_img.get_data())
            if ~(mask_affine == nii_affine).all():
                print("Affine matrices don't match up! " +
                      "Mask index must be adjusted.")
                zero_row = np.zeros(len(mask_index[0]))
                m_idx_4d = np.vstack((mask_index, zero_row))
                shift_matrix = np.linalg.inv(nii_affine).dot(mask_affine)
                adj_m_index = shift_matrix.dot(m_idx_4d).astype(int)[:3]
                mask_index = tuple(map(tuple, adj_m_index))
            masked_vals = nii_vals[mask_index]
            nanmasked_nii = np.ma.masked_array(masked_vals,
                                               np.isnan(masked_vals))
        else:
            raise ValueError("Mask and nifti image are " +
                             "different shapes/sizes.")
        return nanmasked_nii.mean()
    else:
        raise ValueError("Variable mask_image points to non-existent file.")


def str2bool(string):
    """
    Returns boolean value based on string. True strings include 'yes', 'true',
    't', '1', and 'y', along with any of their case variations.
    """
    return string.lower() in ("yes", "true", "t", "1", "y")


def _try_index(list_, val):
    """
    Indexes a list without throwing an error if the value isn't found.
    """
    try:
        return list_.index(val)
    except:
        return False


def determine_spm_contrasts(spm_file):
    mat = scipy.io.loadmat(spm_file, squeeze_me=True, struct_as_record=False)
    spm = mat["SPM"]
    contrast_names = []
    contrast_files = []
    for i, i_con in enumerate(spm.xCon):
        contrast_names.append(str(i_con.name))
        contrast_files.append(str(i_con.Vspm.fname))
        print("%s\t%s" % (str(i), str(i_con.name)))

    sure = False
    while not sure:
        wanted_index = input("Contrast list: ")
        sure = str2bool(input("You sure (yes/no)? "))
    wanted_contrasts = [contrast_names[i] for i in wanted_index]
    return wanted_contrasts


def determine_contrast_files(spm_file, contrast_list):
    mat = scipy.io.loadmat(spm_file, squeeze_me=True, struct_as_record=False)
    spm = mat["SPM"]
    full_file_list = [str(i_con.Vcon.fname) for i_con in spm.xCon]
    full_contrast_list = [str(i_con.name) for i_con in spm.xCon]
    contrast_index = [_try_index(full_contrast_list, contrast) for contrast in
                      contrast_list]
    file_list = [full_file_list[index] for index in contrast_index]
    return file_list


if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise ValueError("Variable nifti_file not given.")
    elif len(sys.argv) == 2:
        nifti_file = sys.argv[1]
        mean = main(nifti_file)
    elif len(sys.argv) == 3:
        nifti_file = sys.argv[1]
        mask_file = sys.argv[2]
        mean = main(nifti_file, mask_file)
    else:
        raise ValueError("Too many inputs given.")
