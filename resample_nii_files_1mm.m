% Resample the lesion masks to 1x1x1 mm.

% Requires the NIfTI and ANALYZE Tools toolbox.

% Created by KD 2022-09-20
% UPD by KD 2023-03-21 to address the change in behavior of 'fileparts',
% which makes the code work in MATLAB R 2022a, but causes error in MATLAB R2020a

% Get ID numbers of VHIS patients (VHIS_ID) in double format

% UPD 2023-03-21: This works in MATLAB R2022a, but not in MATLAB R2020a
% nii_dir = struct2table(dir('nii files'));
% [~, nii_fnames, ~] = fileparts(nii_dir.name(3:length(nii_dir.name)));
% subNum_vec = str2double(nii_fnames)';

% This is a more robust workaround suggested by ChatGPT:
folder_path = [pwd '/nii files'];
file_pattern = '*.nii';
nii_dir = dir(fullfile(folder_path, file_pattern));
[~, nii_fnames, ~] = cellfun(@fileparts, {nii_dir(:).name}, 'UniformOutput', false);

subNum_vec = str2double(nii_fnames);


% Create a folder to store resampled nii images
output_fn = 'nii files resampled_1mm';
if (~exist(output_fn, 'dir')); mkdir(output_fn); end

% Flip each of the lesion masks and store them in the created folder
for iSub = 1:length(subNum_vec)
    old_fname = sprintf('nii files/%04d.nii', subNum_vec(iSub));
    new_fname = sprintf('nii files resampled_1mm/%04d.nii', subNum_vec(iSub));
    reslice_nii(old_fname, new_fname, [1 1 1]);
end