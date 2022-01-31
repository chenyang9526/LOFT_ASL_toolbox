delete(strcat(subject_ID,'_T1w.nii'));
anat_nii = dir(strcat(subject_ID,'_T1w.nii.gz'));
system(strcat(config_info.FSL_path,'/bin/bet',32,anat_nii.name,32,'bet_',anat_nii.name));
system(strcat(config_info.FSL_path,'/bin/fast',32,'bet_',anat_nii.name));
gunzip(strcat('bet_',subject_ID,'_T1w_pve_0.nii.gz'));
gunzip(strcat('bet_',subject_ID,'_T1w_pve_1.nii.gz'));
gunzip(strcat('bet_',subject_ID,'_T1w_pve_2.nii.gz'));
gunzip(strcat(subject_ID,'_T1w.nii.gz'));
