try
    fid = fopen(config_path);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    [config_info,~] = parse_json(str);
    config_info = config_info{1};
catch
    error('Invalid configuration file! The config file must be a json file.');
end
addpath(config_info.NIFTI_path);
addpath(config_info.SPM_path);
if config_info.T1w_flag == 1
    seg_jobfile = {[config_info.CBF_toolbox_path,'/fun_segment_job.m']};
    coreg_jobfile = {[config_info.CBF_toolbox_path,'/fun_coregister_job.m']};
    if config_info.PVC_flag == 1
        if isempty(config_info.FSL_path)
            warning('No FSL_path found in config file. PVC requires FSL./n');
            config_info.FSL_path = input('Enter FSL_path here:\n');
        end
        fsldir = config_info.FSL_path; 
        fsldirmpath = sprintf('%s/etc/matlab',fsldir);
        setenv('FSLDIR', fsldir);
        setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
        path(path, fsldirmpath);
        clear fsldir fsldirmpath;
    end
end