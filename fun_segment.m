% List of open inputs
nrun = 1; % enter the number of runs here
jobs = repmat(seg_jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
clear matlabbatch;
