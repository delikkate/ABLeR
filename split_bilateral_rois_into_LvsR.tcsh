#!/usr/bin/env tcsh

# created by KD
# last update 2023-02-14

# This little piece of code splits each bilateral Brodmann ROI into two hemispheric ROIs using the AFNI function 3dcalc. 

cd ~/Github_code/ABLeR/ROIs_brodmann

foreach i (01 02 03 04 05 06 07 08 09 10 11 17 18 19 20 21 22 23 24 25 26 27 28 29 30 32 34 35 36 37 38 39 40 41 42 43 44 46 46 47 48)
    3dcalc -a roi${i}.nii -prefix roi${i}_LH.nii -expr 'a*step(x)'
    3dcalc -a roi${i}.nii -prefix roi${i}_RH.nii -expr 'a*step(-x)'
end

