<h3>Analysis of Brain Lesions, Revived (ABLeR)</h1>

This repo represents a collection of tools that were developed to characterize brain lesions of patients from the Vietnam Head Injury Study (VHIS) dataset ([Raymont et al. 2011](https://www.frontiersin.org/articles/10.3389/fneur.2011.00015/full)). 
It aims to replace the Analysis of Brain Lesions (ABLe) software ([Solomon et al. 2007](https://www.sciencedirect.com/science/article/pii/S016926070700034X?via%3Dihub); [Makale et al. 2002](https://link.springer.com/article/10.3758/BF03195419)) that is not supported any more.

The `compute_lesion_ROI_overlap.m` script produces the tables containing information about the supplied binary lesion masks, including:
- the overall lesion size (in voxels) and volume (in cm<sup>3</sup>),
- the extent of lesion overlap with different anatomical ROIs (in voxels and cm<sup>3</sup>),
- the percentage of a lesion falling into each anatomical ROI,
- the percentage of each anatomical ROI impacted by a lesion.

This script automates the procedure that can be done manually in [MRIcroGL](https://www.nitrc.org/projects/mricrogl), as demonstrated in [this video tutorial](https://www.youtube.com/watch?v=VfpEv8Y2EP0), and allows to do it for multiple lesion masks in one batch.

The script requires the [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) toolbox for MATLAB to be installed. Download the `spm12` folder and add it to the MATLAB path. If a Mac refuses to run Mex files called by SPM, override the Gatekeeper in the Terminal by typing the following (replace `LOCATION_OF_SPM` with the path where it’s stored):
```
sudo xattr -r -d com.apple.quarantine LOCATION_OF_SPM
sudo find LOCATION_OF_SPM -name \*.mexmaci64 -exec spctl --add {} \;
```

Copy and paste the anatomical atlas (both the `.nii` and the `.xml`) you want to use to `spm12/atlas` or `spm12/tpm`. A few atlases in the MNI space are provided with this code in the subfolder `MNI_atlases`, including:
- three versions of the [AAL atlas](https://search.kg.ebrains.eu/instances/Dataset/f8758eda-483e-45fe-8a88-a1fc806dde18) for SPM:
	- AAL (2&times;2&times;2 mm voxels)
	- AAL3v1 (2&times;2&times;2 mm voxels)
	- AAL3v1_1mm (1&times;1&times;1 mm voxels)
- version of the Brodmann atlas, copied and modified from MRIcroGL’s predecessor [MRIcron](https://www.nitrc.org/projects/mricron). Note that numbers in ROI names in *brodmann.xml* need to be zero-padded, otherwise `spm_atlas.m` will lump all double-digit ROIs into the corresponding single-digit one (so that `roi1`, `roi11`, `roi17`, `roi18` and `roi19` will all be extracted into `roi1`).
To find the atlases supplied with MRIcron or MRIcroGL on a Mac computer, go to Applications &rarr; right-click on the MRIcron/MRIcroGL icon to evoke context menu &rarr; Show Package Contents<br>
&emsp; &rarr; MRIcron: Contents/resources/templates<br>
&emsp; &rarr; MRIcroGL: Contents/Resources/atlas

The script takes three inputs:
1) `atlas_name` &mdash; file name (without the extension) of the anatomical atlas of choice, that needs to be stored in `spm12/atlas` or `spm12/tpm`
2) `lesion folder` &mdash; name of the subfolder containing the binary lesion masks
3) `voxel_volume_cm3` &mdash; e.g., 1-mm isotropic voxel has the volume of (1^3)/10= 0.001 cm<sup>3</sup>; 2-mm isotropic voxel has the volume of (2^3)/1000 = 0.008 cm<sup>3</sup>

Make sure that the anatomical atlas and the lesion masks are in the same space (e.g., MNI) and have the same dimensions and resolution (voxel size). To check voxel size in MRIcroGL, right-click on the atlas name in the Layers window &rarr; Show Header.


\
\
**Supplementary tools:**
- `extract_all_atlas_rois_as_nii.m` &mdash; export anatomical ROIs from an atlas to separate `.nii` files
 - `flip_nii_files_LR.m` &mdash; flip images horizontally, in case the L-R orientation is not recognized correctly by the imaging software/image viewer; uses the `flip_lr.m` function from the [NIfTI and ANALYZE Tools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) toolbox
- `resample_nii_files_1mm.m` &mdash; resample images to 1&times;1&times;1 mm; uses the `reslice_nii.m` function from the [NIfTI and ANALYZE Tools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) toolbox
- `split_bilateral_rois_into_LvsR.sh` &mdash; a piece of T-shell code calling the [AFNI function](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html) `3dcalc` to split the bilateral Brodmann ROIs into two hemispheric ROIs
