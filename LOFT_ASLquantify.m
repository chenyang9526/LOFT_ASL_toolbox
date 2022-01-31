% Usage:
% Step 1: Arrange folder structure to be consistent with that in Doc section 2.8;
% Step 2: Prepare a config json file and copy the path of the config file (see Doc section 2.9);
% Step 3: Prepare a subject json file and store it in the perf folder of each subject folder (see Doc section 2.10);
% Step 4: Add LOFT_ASL_toolbox into your MATLAB paths;
% Step 5: In MATLAB command window, type:
%         LOFT_ASLquantify('your_path_to_config_file/config.json');

function [config_info] = LOFT_ASLquantify(config_path)
if nargin < 1
    error('Missing config file!');
end
%% load config file
fun_load_config;
%% load and unzip OSIPI data
config_info.data_path = [config_info.data_path,'/rawdata'];
cd(config_info.data_path);
dirpath = strcat(config_info.data_path,'/',strtrim(strsplit(config_info.subject_ID,',')));

for i_subject = 1:size(dirpath,2)
    subject_path = string(dirpath(:,i_subject));
    [~,subject_ID,~] = fileparts(subject_path);
    fprintf('Start processing %s \n',subject_ID);
    try
        cd(subject_path);
    catch
        error('Invalid subject folder: %s! Please check your data path in config file and/or your subject folder selection!',string(dirpath(:,i_subject)));
    end
    fun_load
    %% load parameters
    cd perf;
    fun_load_params;
    cd ..;
    %% segment and generate mask based on T1w
    if config_info.T1w_flag == 1
        fprintf('Start segmentation...\n');
        try
            cd anat;
        catch
            error('No T1w data was found. If you do not have T1w data, please set T1w_flag to 0 in config file.');
        end
        try
            evalc('fun_segment');     
            if config_info.PVC_flag == 1
                try
                    evalc('fun_segment_FSL');
                catch
                    error('Partial volume Segmentation failed. Please check your FSL path in config file.');
                end
            end   
        catch
            error('Segmentation failed. Please check your T1w data and/or paths in config file.');
        end
        cd ..;
        fprintf('Segmentation completed\n');
    end
    %% coregister T1w to M0
    if config_info.T1w_flag == 1
        fprintf('Start coregistration of T1w to M0...\n');
        try
            cd anat;
        catch
            error('No anat folder was found. If you do not have T1w data, please set T1w_flag to 0 in config file.');
        end
        try
            evalc('fun_coregister');
            if config_info.PVC_flag == 1
                try
                    evalc('fun_coregister_FSL');
                catch
                    error('Partial volume Segmentation failed. Please check your FSL path in config file.');
                end
            end  
        catch
            error('Coregistration failed. Please check your T1w and M1 data and/or paths in config file.');
        end
        cd ..;
        fprintf('Coregistration completed\n');
    end
    %% expand 4D ASL to 3D
    try
        cd perf;
    catch
        error('No perf folder was found! If you have ASL data please put them in perf folder.');
    end
    try
        expand_nii_scan(dirpath_asl_nii.name);
    catch
        error('Wrong NIFTI_path in config file!');
    end
    cd ..;
    %% MoCo for ASL based on M0
    fprintf('Start motion correction...\n');
    cd perf;
    try
        evalc('fun_moco');
    catch
        error('MoCo failed! Please check your ASL data.');
    end
    cd ..;
    fprintf('Motion correction completed\n');
    %% select files
    cd perf;
    fun_select_files;
    cd ..;
    %% calculate CBF
    fprintf('Start ASL quantification...\n');
    LOFT_CBFquantify(ASL, M0, params);
    fprintf('ASL quantification completed\n');
    %% PVE correction
    fprintf('Start PVE correction...\n');
    if config_info.PVC_flag == 1
        if config_info.T1w_flag == 0
            warning('PVE failed! No T1w data was provided. ');
        else
            fun_PVE_correction;
        end
    end
    fprintf('PVE correction completed\n');
    %% Save files
    fun_save_files;
    %%
    fprintf('Current progress is: %d / %d \n',i_subject,size(dirpath,2));
end
end

