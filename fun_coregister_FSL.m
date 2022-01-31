% List of open inputs
nrun = 1; % enter the number of runs here
coreg_jobfile_FSL = {[config_info.CBF_toolbox_path,'/fun_coregister_job_FSL.m']};
jobs = repmat(coreg_jobfile_FSL, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch;
