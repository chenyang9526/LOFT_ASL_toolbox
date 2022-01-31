
[PATHSTR,Subject_NAME,EXT] = fileparts(string(dirpath(:,i_subject)));
[PATHSTR,NAME,EXT] = fileparts(PATHSTR);
path_subject = strcat(PATHSTR,'/result/',Subject_NAME);
try
    mkdir(path_subject);
catch
    error('Cannot create folder to save result files. Please check your paths.');
end
%%
if config_info.CBF_flag && config_info.mean_flag
    cd perf;
    cbf_nii = dir('*.nii');
    cbf_nii = regexpi({cbf_nii.name},'^AMean_CBF.*','match');
    cbf_nii = [cbf_nii{:}];
    gzip(cbf_nii);
    copyfile(strcat(string(cbf_nii),'.gz'), strcat(path_subject,'/CBF.nii.gz'));
    V = spm_vol(deblank(cbf_nii));
    cbf = spm_read_vols(V{1});
    cd ..;
end

if config_info.CBF_flag && config_info.mean_flag && config_info.PVC_flag && config_info.T1w_flag
    cd perf;
    cbf_gm_nii = dir('*.nii');
    cbf_gm_nii = regexpi({cbf_gm_nii.name},'^PVEc_GM.*','match');
    cbf_gm_nii = [cbf_gm_nii{:}];
    gzip(cbf_gm_nii);
    copyfile(strcat(string(cbf_gm_nii),'.gz'), strcat(path_subject,'/CBF_GMpv.nii.gz'));
    V = spm_vol(deblank(cbf_gm_nii));
    cbf_GMpv = spm_read_vols(V{1});
    
    cbf_wm_nii = dir('*.nii');
    cbf_wm_nii = regexpi({cbf_wm_nii.name},'^PVEc_WM.*','match');
    cbf_wm_nii = [cbf_wm_nii{:}];
    gzip(cbf_wm_nii);
    copyfile(strcat(string(cbf_wm_nii),'.gz'), strcat(path_subject,'/CBF_WMpv.nii.gz'));
    V = spm_vol(deblank(cbf_wm_nii));
    cbf_WMpv = spm_read_vols(V{1});
    cd ..;
end

if config_info.T1w_flag
    cd anat;
    if config_info.PVC_flag
        gm_nii = dir('*.nii');
        gm_nii = regexpi({gm_nii.name},'^rbet_.*_pve_1.nii','match');
        gm_nii = [gm_nii{:}];
        gzip(gm_nii);
        copyfile(strcat(string(gm_nii),'.gz'), strcat(path_subject,'/GM_pv.nii.gz'));
    end
    gm_nii = dir('*.nii');
    gm_nii = regexpi({gm_nii.name},'^rc1.*','match');
    gm_nii = [gm_nii{:}];
    V = spm_vol(deblank(gm_nii));
    mask_gm = spm_read_vols(V{1});
    mask_gm1 = imgaussfilt3(mask_gm,1);
    mask_gm(mask_gm>=0.99) = 1;
    mask_gm(mask_gm~=1) = 0;
    V  = spm_create_vol(V{1});
    V.fname = char(strcat('mask_',string(gm_nii)));
    spm_write_vol(V,mask_gm);
    
    
    if config_info.PVC_flag
        wm_nii = dir('*.nii');
        wm_nii = regexpi({wm_nii.name},'^rbet_.*_pve_2.nii','match');
        wm_nii = [wm_nii{:}];
        gzip(wm_nii);
        copyfile(strcat(string(wm_nii),'.gz'), strcat(path_subject,'/WM_pv.nii.gz'));
    end
    wm_nii = dir('*.nii');
    wm_nii = regexpi({wm_nii.name},'^rc2.*','match');
    wm_nii = [wm_nii{:}];
    V = spm_vol(deblank(wm_nii));
    mask_wm = spm_read_vols(V{1});
    mask_wm(mask_wm>=0.99) = 1;
    mask_wm(mask_wm~=1) = 0;
    mask_wm(mask_wm~=1) = 0;
    SE = strel('cuboid',[2 2 2]);
    mask_wm = imerode(mask_wm,SE);
    V  = spm_create_vol(V{1});
    V.fname = char(strcat('mask_',string(wm_nii)));
    spm_write_vol(V,mask_wm);
    
    gm_nii = dir('*.nii');
    gm_nii = regexpi({gm_nii.name},'^mask_rc1.*','match');
    gm_nii = [gm_nii{:}];
    gzip(gm_nii);
    copyfile(strcat(string(gm_nii),'.gz'), strcat(path_subject,'/GM_mask_lowres.nii.gz'));
    
    wm_nii = dir('*.nii');
    wm_nii = regexpi({wm_nii.name},'^mask_rc2.*','match');
    wm_nii = [wm_nii{:}];
    gzip(wm_nii);
    copyfile(strcat(string(wm_nii),'.gz'), strcat(path_subject,'/WM_mask_lowres.nii.gz'));
    cd ..;
end

cd(path_subject);
if config_info.CBF_flag && config_info.mean_flag && config_info.T1w_flag
    mean_gm = mean(cbf(find(mask_gm == 1)));
    std_gm = std(cbf(find(mask_gm == 1)));
    mean_wm = mean(cbf(find(mask_wm == 1)));
    std_wm = std(cbf(find(mask_wm == 1)));
    
    fileID = fopen('GM_CBF.txt','w');
    fprintf(fileID,'Mean: %.2f\n',mean_gm);
    fprintf(fileID,'Std: %.2f\n',std_gm);
    fclose(fileID);
    
    fileID = fopen('WM_CBF.txt','w');
    fprintf(fileID,'Mean: %.2f\n',mean_wm);
    fprintf(fileID,'Std: %.2f\n',std_wm);
    fclose(fileID);
end

if config_info.CBF_flag && config_info.mean_flag && config_info.T1w_flag && config_info.PVC_flag
    mean_gm_pv = mean(cbf_GMpv(find(mask_gm == 1)));
    std_gm_pv = std(cbf_GMpv(find(mask_gm == 1)));
    mean_wm_pv = mean(cbf_WMpv(find(mask_wm == 1)));
    std_wm_pv = std(cbf_WMpv(find(mask_wm == 1)));
    
    fileID = fopen('GMpv_CBF.txt','w');
    fprintf(fileID,'Mean: %.2f\n',mean_gm_pv);
    fprintf(fileID,'Std: %.2f\n',std_gm_pv);
    fclose(fileID);
    
    fileID = fopen('WMpv_CBF.txt','w');
    fprintf(fileID,'Mean: %.2f\n',mean_wm_pv);
    fprintf(fileID,'Std: %.2f\n',std_wm_pv);
    fclose(fileID);
end
cd ..;