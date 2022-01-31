[subject_ID,~,~] = fileparts(pwd);
[~,subject_ID,~] = fileparts(subject_ID);
anat_nii = dir('*.nii');
anat_path = anat_nii(1).folder;
%anat_nii = regexpi({anat_nii.name},'^[^c\d].*','match');
anat_nii = regexpi({anat_nii.name},char(strcat(subject_ID,'_T1w\.nii$')),'match');
anat_nii = [anat_nii{:}];
source_nii = strcat(anat_path,'/',anat_nii{1},',1');

mask_nii = dir('*.nii');
%mask_nii = regexpi({mask_nii.name},'^c\d.*','match');
mask_nii = regexpi({mask_nii.name},char(strcat('^bet_',subject_ID,'_T1w_pve_\d\.nii$')),'match');
mask_nii = [mask_nii{:}];
mask_nii = strcat(anat_path,'/',mask_nii,',1');
mask_nii = mask_nii.';

cd ..;
cd perf;
m0_nii = dir('*.nii');
perf_path = m0_nii(1).folder;
%m0_nii = regexpi({m0_nii.name},'.*m0scan\.nii','match');
m0_nii = regexpi({m0_nii.name},char(strcat('^',subject_ID,'_m0scan\.nii$')),'match');
m0_nii = [m0_nii{:}];
m0_nii = strcat(perf_path,'/',m0_nii{1},',1');
cd ../anat;

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {m0_nii};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {source_nii};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = mask_nii;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
