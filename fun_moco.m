[subject_ID,~,~] = fileparts(pwd);
[~,subject_ID,~] = fileparts(subject_ID);
ASL = dir('*.nii');
perf_path = ASL(1).folder;
%ASL = regexpi({ASL.name},'^.*asl_\d{4}\.nii$','match');
ASL = regexpi({ASL.name},char(strcat(subject_ID,'_asl_\d{4}\.nii$')),'match');
ASL = [ASL{:}];
ASL = strcat(perf_path,'/',ASL);
ASL = ASL.';

m0_nii = dir('*.nii');
%m0_nii = regexpi({m0_nii.name},'.*m0scan\.nii','match');
m0_nii = regexpi({m0_nii.name},char(strcat('^',subject_ID,'_m0scan\.nii$')),'match');
m0_nii = [m0_nii{:}];
m0_nii = strcat(perf_path,'/',m0_nii);

odd_ASL = [m0_nii;ASL(1:2:end)];
even_ASL = [m0_nii;ASL(2:2:end)];

% % odd images
matlabbatch{1}.spm.spatial.realign.estwrite.data = {odd_ASL};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % 0 reg to first / 1 reg to mean
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch);
clear matlabbatch;

[PATHSTR,NAME,EXT] = fileparts(m0_nii{1});
txt_path = strcat(PATHSTR, '/rp_', NAME,'.txt');
txt_path = cellstr(txt_path);
meanFD.odd = JointMotionParameters(txt_path);
delete(txt_path{1});

% even images
matlabbatch{1}.spm.spatial.realign.estwrite.data = {even_ASL};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % 0 reg to first / 1 reg to mean
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch);
clear matlabbatch;

meanFD.even = JointMotionParameters(txt_path);
delete(txt_path{1});

save('meanFD','meanFD');