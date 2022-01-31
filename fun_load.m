if config_info.T1w_flag == 1
    try
        cd anat;
        dirpath_anat_nii=dir(strcat(subject_ID,'_T1w.nii'));
        if isempty(dirpath_anat_nii)           
            dirpath_anat_nii=dir(strcat(subject_ID,'_T1w.nii.gz'));
            gunzip(dirpath_anat_nii.name);
            dirpath_anat_nii=dir(strcat(subject_ID,'_T1w.nii'));
        end
        cd ..;
    catch
        error('Failed to load T1w data! Check your folder and file names (Should be subject_ID_T1w.nii or .nii.gz). If you do not have T1w data, please set T1w_flag = 0 in config file.');
    end
end

cd perf;
try
    dirpath_asl_json=dir(strcat(subject_ID,'_asl.json'));
catch
    error('json file for perf data is not provided!');
end
try
    dirpath_asl_nii=dir(strcat(subject_ID,'_asl.nii'));
    if isempty(dirpath_asl_nii)
        dirpath_asl_nii=dir(strcat(subject_ID,'_asl.nii.gz'));
        gunzip(dirpath_asl_nii.name);
        dirpath_asl_nii=dir(strcat(subject_ID,'_asl.nii'));
    end
catch
    error('Failed to load ASL data! Check your folder and file names (Should be subject_ID_asl.nii or .nii.gz).');
end
if size(dirpath_asl_nii,1) > 1
    error('Multiple ASL data was found! Check your ASL data. ASL data should be 4D.');
end
try
    dirpath_m0_nii=dir(strcat(subject_ID,'_m0scan.nii'));
    if isempty(dirpath_m0_nii) 
        dirpath_m0_nii=dir(strcat(subject_ID,'_m0scan.nii.gz'));
        gunzip(dirpath_m0_nii.name);
        dirpath_m0_nii=dir(strcat(subject_ID,'_m0scan.nii'));
    end
catch
    error('Failed to load M0 data! Check your folder structure.');
end
if size(dirpath_asl_nii,1) > 1
    error('Multiple M0 data was found! Check your M0 data (Should be subject_ID_m0scan.nii or .nii.gz).');
end
cd ..