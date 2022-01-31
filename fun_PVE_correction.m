
cd anat;
masks = dir('*.nii');
GM_mask = regexpi({masks.name},'^rbet.*pve_1\.nii$$','match');
GM_mask = [GM_mask{:}];
if isempty(GM_mask)
    warning('No GM mask file were found! Please check if you have a file named rc1...nii in anat folder!');
end
GM_mask = spm_vol(deblank(GM_mask));
GM_mask = spm_read_vols(GM_mask{1});
GM_mask(isnan(GM_mask)) = 0;
GM_mask(isinf(GM_mask)) = 0;
GM_mask(find(GM_mask>1)) = 1;
GM_mask(find(GM_mask<0)) = 0;

WM_mask = regexpi({masks.name},'^rbet.*pve_2\.nii$$','match');
WM_mask = [WM_mask{:}];
if isempty(WM_mask)
    warning('No GM mask file were found! Please check if you have a file named rc1...nii in anat folder!');
end
WM_mask = spm_vol(deblank(WM_mask));
WM_mask = spm_read_vols(WM_mask{1});
WM_mask(isnan(WM_mask)) = 0;
WM_mask(isinf(WM_mask)) = 0;
WM_mask(find(WM_mask>1)) = 1;
WM_mask(find(WM_mask<0)) = 0;

CSF_mask = regexpi({masks.name},'^rbet.*pve_0\.nii$$','match');
CSF_mask = [CSF_mask{:}];
if isempty(CSF_mask)
    warning('No GM mask file were found! Please check if you have a file named rc1...nii in anat folder!');
end
CSF_mask = spm_vol(deblank(CSF_mask));
CSF_mask = spm_read_vols(CSF_mask{1});
CSF_mask(isnan(CSF_mask)) = 0;
CSF_mask(isinf(CSF_mask)) = 0;
CSF_mask(find(CSF_mask>1)) = 1;
CSF_mask(find(CSF_mask<0)) = 0;
cd ..

cd perf;
mCBF = dir('*.nii');
mCBF = regexpi({mCBF.name},'^AMean_CBF.*\.nii$$','match');
mCBF = [mCBF{:}];
if isempty(mCBF)
    error('No mean CBF file were found! Please check if you have a file named AMean_CBF...nii in perf folder!');
end
VCBF = spm_vol(deblank(mCBF));
[pth,nm,xt] = fileparts(deblank(VCBF{1}.fname));
mCBF = spm_read_vols(VCBF{1});
mCBF(isnan(mCBF)) = 0;
mCBF(isinf(mCBF)) = 0;
%%
mCBF_GM_PVEc = zeros(size(mCBF));
mCBF_WM_PVEc = zeros(size(mCBF));
for z = 1:size(mCBF,3)
    for x = 1:size(mCBF,1)
        for y = 1:size(mCBF,2)
            if WM_mask(x,y,z) == 0 && GM_mask(x,y,z) == 0
                continue;
            end
            if GM_mask(x,y,z) >= 0.999 || (WM_mask(x,y,z) == 0 && GM_mask(x,y,z) ~= 0)
                mCBF_GM_PVEc(x,y,z) = mCBF(x,y,z);
                mCBF_WM_PVEc(x,y,z) = 0;
                continue;
            elseif WM_mask(x,y,z) >= 0.999 || (GM_mask(x,y,z) == 0 && WM_mask(x,y,z) ~= 0)
                mCBF_WM_PVEc(x,y,z) = mCBF(x,y,z);
                mCBF_GM_PVEc(x,y,z) = 0;
                continue;
            end
            pos_x = [x-1,x-1,x-1,x  ,x  ,x  ,x+1,x+1,x+1 ];
            pos_y = [y-1,y  ,y+1,y-1,y  ,y+1,y-1,y  ,y+1 ];
            % check for invalid pos (exceed boundaries)
            temp_x = find(pos_x < 2 | pos_x > size(mCBF,1)-1);
            temp_y = find(pos_y < 2 | pos_y > size(mCBF,2)-1);
            pos_ind_to_clear = union(temp_x,temp_y);
            pos_x(pos_ind_to_clear) = [];pos_y(pos_ind_to_clear) = [];
            % check for invalid pos (not within brain mask)
            pos_ind_to_clear1 = [];
            pos_ind_to_clear2 = [];
            if length(pos_x) < 2
                continue;
            else
                n = length(pos_x);
                PV_matrix = zeros(n,2);
                CBF_vector = zeros(n,1);
                for i_n = 1:n
                    GM_PV = GM_mask(pos_x(i_n),pos_y(i_n),z);
                    WM_PV = WM_mask(pos_x(i_n),pos_y(i_n),z);
                    CSF_PV = CSF_mask(pos_x(i_n),pos_y(i_n),z);
                    PV_sum = GM_PV+WM_PV+CSF_PV;
                    if PV_sum == 0
                        continue;
                    end
                    GM_PV = GM_PV/PV_sum;
                    WM_PV = WM_PV/PV_sum;
                    PV_matrix(i_n,1) = GM_PV;
                    PV_matrix(i_n,2) = WM_PV;
                    CBF_vector(i_n) = mCBF(pos_x(i_n),pos_y(i_n),z);
                end
                if length(find(CBF_vector > 0)) < 2
                    continue;
                end
                [U,S,V] = svd(PV_matrix);
                pinv_PV_matrix = pinv(PV_matrix,max(S(:))*0.15);
                CBF_PVEc = pinv_PV_matrix*CBF_vector;
                if CBF_PVEc(1) >= 0.1 && CBF_PVEc(1) < CBF_PVEc(2)
                    mCBF_GM_PVEc(x,y,z) = CBF_PVEc(1);
                    mCBF_WM_PVEc(x,y,z) = 0;
                else
                    mCBF_GM_PVEc(x,y,z) = CBF_PVEc(1);
                    mCBF_WM_PVEc(x,y,z) = CBF_PVEc(2);
                end
                mCBF_GM_PVEc(x,y,z) = CBF_PVEc(1);
                mCBF_WM_PVEc(x,y,z) = CBF_PVEc(2);
            end
        end
    end
end
GM_mask_temp = zeros(size(GM_mask));
GM_mask_temp(find(GM_mask)<0.05) = 1;
WM_mask_temp = zeros(size(WM_mask));
WM_mask_temp(find(WM_mask)<0.05) = 1;
mCBF_GM_PVEc = mCBF_GM_PVEc.*(1-GM_mask_temp);
mCBF_WM_PVEc = mCBF_WM_PVEc.*(1-WM_mask_temp);
mCBF_GM_PVEc(find(mCBF_GM_PVEc)<0) = 0;
mCBF_WM_PVEc(find(mCBF_WM_PVEc)<0) = 0;

VGM=VCBF{1};
VGM.fname    = fullfile(pth,['PVEc_GM_' nm xt]);
VGM.dt = [16 0];
VGM  = spm_create_vol(VGM);
VGM = spm_write_vol(VGM,mCBF_GM_PVEc);

VWM=VCBF{1};
VWM.fname    = fullfile(pth,['PVEc_WM_' nm xt]);
VWM.dt = [16 0];
VWM  = spm_create_vol(VWM);
VWM = spm_write_vol(VWM,mCBF_WM_PVEc);
cd ..