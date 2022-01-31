
M0 = dir('*.nii');
%M0 = regexpi({M0.name},'.*m0scan\.nii$','match');
M0 = regexpi({M0.name},char(strcat(subject_ID,'_m0scan\.nii$')),'match');
M0 = [M0{:}];
if isempty(M0)
    error('No M0 files were found! Your M0 files must have suffix of m0scan.nii.gz or m0scan.nii!');
end
[~,M0_idx] = min(cellfun('length', M0));
if length(M0_idx)>1
    error('More than 1 M0 files were selected. Please delete the duplicated M0.');
end
M0={['perf/',M0{M0_idx}]};

ASL = dir('*.nii');
%ASL = regexpi({ASL.name},'^[^r].*asl_\d{4}\.nii$','match');
ASL = regexpi({ASL.name},char(strcat('^',subject_ID,'_asl_\d{4}\.nii$')),'match');
ASL = [ASL{:}];
if isempty(ASL)
    error('No ASL files were found! Your ASL files must have suffix of asl.nii.gz or asl.nii!');
end
ASL = strcat('perf/',ASL);
