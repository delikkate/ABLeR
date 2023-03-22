%% Create separate .nii masks for all anatomical atlas ROIs
% Based on this script: https://www.nemotos.net/?p=1083

% Created by KD 2022-12-08
% Last update by KD 2023-03-21: put created ROIs in a new folder

% Select the anatomical atlas
atlas_name = 'brodmann'; %'AAL3v1'

% Load the anatomical atlas using SPM
atlas = spm_atlas('load', atlas_name);
ROI_names = struct2table(atlas.labels).name';


% Create a folder for new .nii files
output_fn = sprintf('ROIs_%s', atlas_name);
if (~exist(output_fn, 'dir')); mkdir(output_fn); end

for iROI = 1:length(atlas.labels)
    ROI_mask = spm_atlas('mask', atlas, atlas.labels(iROI).name); % make sure to zero-pad roi numbers in /spm12/atlas/brodmann.xml, otherwise SPM lumps all double-digit rois into the corresponding single-digit roi
    ROI_mask.fname = sprintf('%s/%s.nii', output_fn, ROI_names{iROI});
    spm_write_vol(ROI_mask, spm_read_vols(ROI_mask));
end

