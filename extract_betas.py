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


def main(nifti_image, mask_image=""):
    if os.path.isfile(nifti_image):
        nii_img = nb.load(nifti_image)
        nii_vals = np.array(nii_img.get_data())
        nii_header = nii_img.get_header()
        nii_shape = nii_header.get_data_shape()
        nii_affine = nii_header.get_qform()
    else:
        raise ValueError("Variable nifti_image points to non-existent file.")

    if os.path.isfile(mask_image):
        mask_img = nb.load(mask_image)
        mask_header = mask_img.get_header()
        mask_shape = mask_header.get_data_shape()
        mask_affine = mask_header.get_qform()

        if mask_shape == nii_shape:
            mask_index = np.where(mask_img.get_data())
            if (mask_affine == nii_affine).all():
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
    elif mask_image:
        raise ValueError("Variable mask_image points to non-existent file.")
    else:
        nanmasked_nii = np.ma.masked_array(nii_vals, np.isnan(nii_vals))
        return nanmasked_nii.mean()

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
