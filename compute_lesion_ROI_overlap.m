% This script produces the tables containing information about the supplied
% binary lesion masks, including:
% 	- the overall lesion size (in voxels) and volume (in cm3),
%   - the extent of lesion overlap with different anatomical ROIs (in
%   voxels and cm3),
%   - the percentage of a lesion falling into each anatomical ROI,
%   - the percentage of each anatomical ROI impacted by a lesion.
%
% The script requires the SPM12 toolbox for MATLAB to be installed. Copy
% and paste the anatomical atlas you want to use to spm12/atlas or
% spm12/tpm.

% Created by KD
% Update log:
% 2022-12-08: fixed issue with extracting single-digit ROIs from Brodmann atlas;
%             added saving workspace into .mat file
% 2023-03-21: introduced a workaround robust to the change in 'fileparts'
%             behavior that shouldn't throw error in older MATLAB versions
%             (tested with MATLAB R2020a and R2022a)


%% Define your variables
% Load an anatomical atlas from SPM (should be stored in spm12/atlas or spm12/tpm)
% Make sure that the atlas has the same dimensions and voxel size as the lesion masks!
atlas_name = 'brodmann'; % 'AAL' | 'AAL3v1' | 'AAL3v1_1mm'

% Indicate path to lesion masks
lesion_folder = 'nii files resampled_1mm'; % 'nii files'

% Enter the volume of a single voxel (to find its dimensions in MRIcroGL: right-click on the image name from the "Layers" list --> Spacing)

% - For Brodmann and AAL3v1_1mm atlases:
voxel_volume_cm3 = (1/10)^3;    % = (1^3)/1000 = 0.001 cm3 = 1 mm3

% - For AAL v.1:
%voxel_volume_cm3 = (2/10)^3;     % = (2^3)/1000 = 0.008 cm3 = 8 mm3


% Start the timer to see how long it takes the output table to be printed
tStart = tic;


%% Load lesion and atlas data
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

% Load the anatomical atlas using SPM
atlas = spm_atlas('load', atlas_name);
ROI_names = struct2table(atlas.labels).name';


%% Preallocate cells arrays for output tables
output_atlasROIsize = cell(length(atlas.labels)+1, 4);
output_atlasROIsize(1,:) = [{'ROI_label'}, {'ROI_index'}, {'Size_in_voxels'}, {'Volume_cm3'}];

output_lesionSize = cell(length(nii_fnames), 5);
output_lesionSize(1,:) = [{'VHIS_ID'}, {'Size_in_voxels'}, {'Volume_cm3'}, {'Voxels_not_covered_by_atlas'}, {'Percentage_of_lesion_not_covered_by_atlas'}];

output_LesionROIvoxelNum = cell(length(nii_fnames)+1, length(atlas.labels)+1);  % number of lesioned voxels per ROI for each subject

output_LesionROIcm3 = cell(length(nii_fnames)+1, length(atlas.labels)+1); % lesion volume per ROI for each subject

output_PercentROIImpactedByLesion = cell(length(nii_fnames)+1, length(atlas.labels)+1);   % percent of ROI impacted by lesion

output_PercentLesionInROI = cell(length(nii_fnames)+1, length(atlas.labels)+1);   % percent of lesion falling into a given ROI




%% Loop through all individual lesion maps (VHIS_ID) and all atlas ROIs
for iSub = 1:length(subNum_vec)

    %% Step 1: Load the binary .nii lesion mask and find the indices of non-zero voxels
    lesion_fname = [lesion_folder sprintf('/%04d.nii', subNum_vec(iSub))];
    lesion.hdr = spm_vol(lesion_fname);
    lesion.dat = spm_read_vols(lesion.hdr);
        
    lesion_ones = find(lesion.dat == 1);   % find indices of voxels included in the lesion mask
    
    output_lesionSize{iSub+1,1} = subNum_vec(iSub);
    output_lesionSize{iSub+1,2} = length(lesion_ones);  % compute lesion size in voxels
    output_lesionSize{iSub+1,3} = round((output_lesionSize{iSub+1,2} * voxel_volume_cm3), 2);  % compute lesion volume in cm3
    

    for iROI = 1:length(atlas.labels)  % loop through all atlas ROIs
        
        output_atlasROIsize{iROI+1,1} = atlas.labels(iROI).name;   % print anatomical label of atlas ROI
        output_atlasROIsize{iROI+1,2} = atlas.labels(iROI).index;   % print numeric index of atlas ROI
        
        % Print ROI name in the first row
        output_LesionROIvoxelNum{1, iROI+1} = struct2table(atlas.labels).name{iROI};
        output_LesionROIcm3{1, iROI+1} = struct2table(atlas.labels).name{iROI};
        output_PercentROIImpactedByLesion{1, iROI+1} = struct2table(atlas.labels).name{iROI};
        output_PercentLesionInROI{1, iROI+1} = struct2table(atlas.labels).name{iROI};
        

        %% Step 2: Load the binary ROI mask and find the indices of non-zero voxels
        ROI_mask = spm_atlas('mask', atlas, atlas.labels(iROI).name);
        ROI_ones = find(ROI_mask.dat == 1); % find indices of voxels that are included in the ROI mask
        output_atlasROIsize{iROI+1,3} = length(ROI_ones); % find how many voxels are in the ROI mask
        output_atlasROIsize{iROI+1,4} = round((output_atlasROIsize{iROI+1,3} * voxel_volume_cm3), 2); % find the volume of the ROI mask in cm3

     
        %% Step 3: Find the overlapping non-zero voxels in the ROI mask (ROI_mask.dat) and the lesion mask (dat)
        overlapping_voxels = intersect(ROI_ones, lesion_ones); % find indices of overlapping voxels (included both in the ROI mask and the lesion mask)
        output_LesionROIvoxelNum{iSub+1, iROI+1} = length(overlapping_voxels);
        output_LesionROIcm3{iSub+1, iROI+1} = round(output_LesionROIvoxelNum{iSub+1, iROI+1} * voxel_volume_cm3, 2);
        output_PercentROIImpactedByLesion{iSub+1, iROI+1} = round((length(overlapping_voxels) / length(ROI_ones) * 100), 2);
        output_PercentLesionInROI{iSub+1, iROI+1} = round((length(overlapping_voxels) / length(lesion_ones) * 100), 2);

        
        % Print VHIS_ID in the first column       
        output_LesionROIvoxelNum{1, 1} = 'VHIS_ID';
        output_LesionROIvoxelNum{iSub+1, 1} = subNum_vec(iSub);
        %outputData{iSub, 1} = sprintf('%04d', subNum_vec(iSub));

        output_LesionROIcm3{1, 1} = 'VHIS_ID';
        output_LesionROIcm3{iSub+1, 1} = subNum_vec(iSub);

        output_PercentROIImpactedByLesion{1, 1} = 'VHIS_ID';
        output_PercentROIImpactedByLesion{iSub+1, 1} = subNum_vec(iSub);

        output_PercentLesionInROI{1, 1} = 'VHIS_ID';
        output_PercentLesionInROI{iSub+1, 1} = subNum_vec(iSub);        


    end

    output_lesionSize{iSub+1,4} = output_lesionSize{iSub+1,2} - ...
        sum(cell2mat(output_LesionROIvoxelNum(iSub+1, 2:end))); % 'Voxels_not_covered_by_atlas'

    output_lesionSize{iSub+1,5} = round((output_lesionSize{iSub+1,4} / output_lesionSize{iSub+1,2} * 100), 2); % 'Percentage_of_lesion_not_covered_by_atlas'

end

%% Step 4. Print the output tables

%xlswrite('VHIS_ALL_sizeInVoxels.xlsx', outputData); % doesn't work on a Mac
writecell(output_LesionROIvoxelNum, sprintf('%s_VHIS_numberLesionedVoxels.xlsx', atlas_name));
writecell(output_LesionROIcm3, sprintf('%s_VHIS_lesionVolume_cm3.xlsx', atlas_name));
writecell(output_atlasROIsize, sprintf('%s_AnatomicalROISize.xlsx', atlas_name));
writecell(output_lesionSize, sprintf('%s_VHIS_LesionSize_nii_files.xlsx', atlas_name));
writecell(output_PercentROIImpactedByLesion, sprintf('%s_VHIS_Percentage_of_ROI_Impacted_By_Lesion.xlsx', atlas_name));
writecell(output_PercentLesionInROI, sprintf('%s_VHIS_Percentage_of_Lesion_In_ROI.xlsx', atlas_name));


% Stop the timer
tEnd_sec = toc(tStart);
tEnd_min = tEnd_sec/60; % the script takes 8.3 min to run on the office MBP


%% Step 5. Save workspace variables for debugging
save(atlas_name)


%% Step 6. Quality control
% Make sure the generated summary table produces the same result as the
% ones generated manually with MRIcroGL. Visually inspect the overlap
% between lesion masks and anatomical ROIs as a sanity check.





