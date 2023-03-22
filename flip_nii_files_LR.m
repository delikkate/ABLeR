% Change the left-right image orientation of lesion masks for them to be
% properly loaded in MRIcroGL. This is a quicker fix than finding what flag
% to set in the header to read the image correctly. Note that in order to
% work with SPM we will use the original images from the "nii files"
% folder, since SPM identifies image orientation correctly.

% Requires the NIfTI and ANALYZE Tools toolbox.

% Created by KD 2022-09-14
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


% Create a folder to store flipped nii images
output_fn = 'nii files flipped';
if (~exist(output_fn, 'dir')); mkdir(output_fn); end

% Flip each of the lesion masks and store them in the created folder
for iSub = 1:length(subNum_vec)
    old_fname = sprintf('nii files/%04d.nii', subNum_vec(iSub));
    new_fname = sprintf('nii files flipped/%04d.nii', subNum_vec(iSub));
    flip_lr(old_fname, new_fname);
end
