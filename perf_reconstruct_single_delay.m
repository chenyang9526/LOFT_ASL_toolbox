% function [gcbf,cbfdat] = perf_reconstruct_2D((Filename, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, ThreshFlag, threshold, CBFFlag, MeanFlag, Delaytime, Labeltime, Slicetime, T1b, PASLMo, optionPCA, ASLscaling, M0scaling, diffs)
%
% This MATLAB function is to reconstruct the raw perfusion images from EPI images by the subtraction between labelled images
% and control images. Quantified CBF images can also be reconstructed by select the option. It is based on SPM2 and
% MATLAB(5.3 or above) on Redhat Linux 9 or WindowsXP/2000.
%
% The MATLAB code of this function will be found in: http://cfn.upenn.edu/perfusion/software.htm
%
% All the images are 3D SPM ANALYZE formatted (.img and .hdr). All the results are also saved in SPM ANALYZE format;
% The labelled and control images should be the data after motion correction.
%
% The method used here are based on the "simple subtraction", "surround subtraction" and "sinc subtraction" approaches described in
% Aguirre GK et al (2002) Experimental design and the relative sensitivity of perfusion and BOLD fMRI, NeuroImage. 15:488-500.
%
% BOLD data (or whatever the underlying pulse sequence that was used) are generated in addition to the perfusion data
%
% for CASL,
% CBF data are calculated according to the formula from
% Wang J, Alsop DC, et al. (2003) Arterial transit time imaging with flow encoding arterial spin tagging (FEAST).
% Magn Reson Med. 50:599-607. Page600, formula [1]
% CBF_CASL (ml/100g/min) = 60*100*M*¦Ë*R/(2*alp*Mo*(exp(-w*R)-exp(-(t+w)*R))
% where M = raw ASL signal, ¦Ë = blood/tissue water partition coefficient, R =longitudinal relaxation rate of blood,
%       alp = tagging efficiency, Mo =  equilibrium magnetization of brain,
%       w = post-labeling delay, t = duration of the labeling pulse,
% and we use the assumed parameters for calculation as ¦Ë=0.9g/ml,
% for 3T, alp=0.68, T1b=1490ms, R=1/T1b=0.67sec-1.
% for 1.5T, alp=0.71, T1b=1200ms, R=1/T1b=0.83sec-1.
%
% for PASL,
% CBF data are calculated according to the fopwdrmula from
% Wang J, Aguirre GK, et al. (2003) Arterial Spin Labeling Perfusion fMRI With Very Low Task Frequency
% Magn Reson Med. 49:796-802. Page798, formula [1]
% CBF_PASL (ml/100g/min) = 60*100*M*¦Ë/(2*alp*Mo*t*exp(-(t+w)*R))
% where M = raw ASL signal, ¦Ë = blood/tissue water partition coefficient, R =longitudinal relaxation rate of blood,
%       alp = tagging efficiency, Mo =  equilibrium magnetization of brain,
%       w = post-labeling delay, t = duration of the labeling pulse,
% and we use the assumed parameters for calculation as ¦Ë=0.9g/ml,
% for 3T, alp=0.95, T1b=1490ms, R=1/T1b=0.67sec-1.
% for 1.5T, alp=0.95, T1b=1200ms, R=1/T1b=0.83sec-1.
%
%  Inputs:
%    Firstimage - integer variable indicating the type of first image
%    - 0:control; 1:labeled
%   Select raw images (*.img, images in a order of control1.img, label1.img, control2.img, label2.img,....;
%   or images in a order of label1.img, control1.img, label2.img, control2.img, .... )
%
%    SubtractionType - integer variable indicating which subtraction method will be used
%    -0: simple subtraction; 1: surround subtraction;2: sinc subtractioin.
%    for CASL, suppose Perfusion = Control - Label;
%    if the raw images is: (C1, L1, C2, L2, C3...)
%     the simple subtraction is: (C1-L1, C2-L2...)
%     the surround subtraction is: ((C1+C2)/2-L1, (C2+C3)/2-L2,...)
%     the sinc subtraction is: (C3/2-L1, C5/2-L2...)
%    if the raw images is: (L1, C1, L2, C2...)
%     the simple subtraction is: (C1-L1, C2-L2...)
%     the surround subtraction is: (C1-(L1+L2)/2, C2-(L2+L3)/2,...)
%     the sinc subtraction is: (C1-L3/2, C2-L5/2...)
%
%    for PASL, suppose Perfusion = Label - Control;
%    if the raw images is: (C1, L1, C2, L2, C3...)
%     the simple subtraction is: (L1-C1, L2-C2...)
%     the surround subtraction is: (L1-(C1+C2)/2, L2-(C2+C3)/2,...)
%     the sinc subtraction is: (L1-C3/2, L2-C5/2...)
%    if the raw images is: (L1, C1, L2, C2...)
%     the simple subtraction is: (L1-C1, L2-C2...)
%     the surround subtraction is: ((L1+L2)/2-C1, (L2+L3)/2-C2,...)
%     the sinc subtraction is: (L3/2-C1, L5/2-C2...)
%
%    ThreshFlag - integer variable indicating whether Threshold EPI image series
%    - 0:no Threshold; 1:Threshold
%
%    CBFFlag - integer variable indicating whether CBF images are produced
%    - 0:no CBF images; 1: produced CBF images
%
%    MeanFlag - integer variable indicating whether mean image of all perfusion images are produced
%    - 0:no mean image; 1: produced mean image
%
%
%  Outputs:
%    BOld Images: Bold_*.img,Bold_*.hdr;  Mean_Bold.img, Mean_Bold.hdr;
%    Perfusion Images: Perf_*.img, Perf_*.hdr; Mean_Perf.img, Mean_Perf.hdr;
%    CBF Images: CBF_*.img, CBF_*.hdr; Mean_CBF.img, Mean_CBF.hdr;
%
%  By H.Y. Rao & J.J. Wang, @CFN, UPenn Med. 07/2004.
%
%
%  Modified by A. Federspiel Uni Bern, Switzerland 07/2008.
%
%
%  adapted for spm8 by A. Federspiel Uni Bern, Switzerland 09/2009.  add ASLscaling by Danny 08/2017
%
function [gcbf,cbfdat, mCBF, Mask] = perf_reconstruct_single_delay(Filename, FieldStrength, ASLType, BSScaling, FirstimageType,...
    SubtractionOrder, SubtractionType, ThreshFlag, threshold, FrameFlag,PerfusionFlag, CBFFlag, MeanFlag, BOLDFlag, AnatFlag, Delaytime, Labeltime, Slicetime, PASLMo, optionPCA, ASLscaling, M0scaling, diffs)

%fprintf('...starting CBF computing...\n');  tic;


if SubtractionType==2,
    Timeshift = 0.5;
end;


if AnatFlag==1,
    cd anat;
    mask_nii = dir('*.nii');
    mask_nii = regexpi({mask_nii.name},'^rc\d.*','match');
    mask_nii = [mask_nii{:}];
    mask_nii = char(mask_nii.');
    if length(mask_nii) == 0
        mask_nii = spm_select(Inf,'image','Select mask images');
    end
    cd ..;
end


if FieldStrength==0,
    R = 0.7407; % T1b = 1350ms at 1.5T based on white paper 
elseif FieldStrength==1,
    R = 0.6061;  % 1/1.650 at 3T based on white paper
elseif FieldStrength==2,
    R = 0.4792; % 1/2.087 at 7T based on Zhang, Xiaofei, et al. "In vivo blood T1 measurements at 1.5T, 3T, and 7T." MRM 2013
end;

if ASLType==1, %PASL
    alp = 0.98;   %tagging efficiency
end;

if ASLType==0, %CASL
    if FieldStrength==0, alp=0.71; end;
    if FieldStrength==1, alp=0.68; end;
    if FieldStrength==2, alp=0.68; end;
end;
if ASLType==2 || ASLType==3, %pCASL, 1.5T and 3T were based on white paper
    if FieldStrength==0, alp=0.85; end;
    if FieldStrength==1, alp=0.85; end; 
    if FieldStrength==2, alp=0.7; end; %need to confirm for 7T
end;
alp = alp*BSScaling;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the main program
% Map images
tmp=Filename(1,:);
prefix=[tmp(1:length(tmp)-3) 'txt'];
clear tmp;

V=spm_vol(deblank(Filename));


%if ASLType==1 ,
VMo = spm_vol(deblank(PASLMo));
PASLModat = zeros([VMo.dim(1:2) 1]);
%end;

perfnum=fix(length(V)/2);



if length(V)==0, error('no raw img files was selected'); end;
if rem(length(V),2)==1, warning('the number of raw img files is not even, last img is ignored'); end;


% Create output images...
VCBF=V(1:perfnum);
VMCBF = V(1);

if BOLDFlag == 1
    VB=V(1:perfnum);
    VMC = V(1);
end

VO = V(1:perfnum);
VMP = V(1);
for k=1:length(V),
    [pth,nm,xt] = fileparts(deblank(V(k).fname));
    if SubtractionType==0,
        VO(k).fname = fullfile(pth,['Perf_0' nm xt]);
        if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_0_' nm xt]);end;
    end;
    if SubtractionType==1,
        VO(k).fname = fullfile(pth,['Perf_1' nm xt]);
        if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_1_' nm xt]);end;
    end;
    if SubtractionType==2,
        VO(k).fname = fullfile(pth,['Perf_2' nm xt]);
        if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_2_' nm xt]);end;
    end;
    
    if BOLDFlag == 1
        VB(k).fname    = fullfile(pth,['Bold_' nm xt]);
        VB(k).dt = [16 0];
    end
    VO(k).dt = [16 0];
    
    if CBFFlag == 1
        VCBF(k).dt = [16 0];
    end
end;


dat  = zeros([VO(1).dim(1:3) length(V)]);

temp=load_untouch_nii(Filename(1,:));
dat = spm_read_vols(V);
dat=dat./temp.hdr.dime.scl_slope;


% read the Mo data for PASL
tempMo=load_untouch_nii(PASLMo);
PASLModat = spm_read_vols(VMo);
PASLModat=PASLModat./tempMo.hdr.dime.scl_slope;
PASLModat=double(PASLModat)./M0scaling;



if AnatFlag == 1
    cd anat
    Mask = spm_vol(deblank(mask_nii));
    Mask = spm_read_vols(Mask);
    Mask = sum(Mask(:,:,:,1:3),4);
    Mask(find(Mask >= 0.5)) = 1;
    Mask(find(Mask < 0.5)) = 0;
    Mask = imfill(Mask,'holes');
    cd ..
elseif ThreshFlag ==1
    Mask=zeros(size(PASLModat));
    Mask(find(PASLModat>(threshold*(max(PASLModat(:))))))=1;
else
    Mask=ones(size(PASLModat));
end

for k=1:length(V),
    dat(:,:,:,k) =(double(dat(:,:,:,k)).*Mask)./ASLscaling;
end;

if FrameFlag == 1
    for k=1:perfnum,
        if BOLDFlag == 1
            VB(k)  = spm_create_vol(VB(k));
        end
        if PerfusionFlag == 1
            VO(k)  = spm_create_vol(VO(k));
        end
        if CBFFlag == 1
            VCBF(k)  = spm_create_vol(VCBF(k));
        end
    end;
end

if BOLDFlag == 1 && FrameFlag == 1
    bdat = zeros([VB(1).dim(1:3) perfnum]);
    for k=1:perfnum,
        bdat(:,:,:,k) = (dat(:,:,:,2*k-1) + dat(:,:,:,2*k))/2;
    end;
    % Write Bold and perfusion image...
    for k=1:perfnum,
        VB(k) = spm_write_vol(VB(k),bdat(:,:,:,k));
        spm_progress_bar('Set',k);
    end;
    clear bdat;
end

cdat = zeros([VO(1).dim(1:3) perfnum]);
ldat = zeros([VO(1).dim(1:3) perfnum]);

% define the control and label images...
for k=1:length(V),
    if SubtractionOrder==0,
        if rem(k,2)== 1, ldat(:,:,:,(k+1)/2) = dat(:,:,:,k); end;
        if rem(k,2)== 0, cdat(:,:,:,k/2) = dat(:,:,:,k); end;
    end;
    if SubtractionOrder==1,
        if rem(k,2)== 1, cdat(:,:,:,(k+1)/2) = dat(:,:,:,k); end;
        if rem(k,2)== 0, ldat(:,:,:,k/2) = dat(:,:,:,k); end;
    end;
end;

clear dat;



pdat = zeros([VO(1).dim(1:3) perfnum]);

if optionPCA==0
    % do the simple subtraction...
    if SubtractionType==0,
        for k=1:perfnum,
            pdat(:,:,:,k) = cdat(:,:,:,k) - ldat(:,:,:,k);
            spm_progress_bar('Set',k);
        end;
    end;
    
    %[perfnum mean(pdat(:))]
    
    % do the linear interpolation...
    if SubtractionType==1,
        %     linear_cdat=zeros([VB(1).dim(1:3) 2*perfnum]);
        %     linear_ldat=zeros([VB(1).dim(1:3) 2*perfnum]);
        linear_cdat=zeros([VO(1).dim(1:3) 2*perfnum]);
        linear_ldat=zeros([VO(1).dim(1:3) 2*perfnum]);
        pnum=1:perfnum;
        lnum=1:0.5:perfnum;
        for x=1:V(1).dim(1),
            for y=1:V(1).dim(2),
                for z=1:V(1).dim(3),
                    cdata = zeros(1,perfnum);
                    ldata = zeros(1,perfnum);
                    linear_cdata = zeros(1,length(V));
                    linear_ldata = zeros(1,length(V));
                    for k=1:perfnum,
                        cdata(k) = cdat(x,y,z,k);
                        ldata(k) = ldat(x,y,z,k);
                    end;
                    linear_cdata=interp1(pnum,cdata,lnum);
                    linear_ldata=interp1(pnum,ldata,lnum);
                    for k=1:2*perfnum-1,
                        linear_cdat(x,y,z,k)= linear_cdata(k);
                        linear_ldat(x,y,z,k)= linear_ldata(k);
                    end;
                end;
            end;
        end;
        
        
        % do the surround subtraction....
        if FirstimageType ==1;
            pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
            spm_progress_bar('Set',1);
            for k=2:perfnum,
                pdat(:,:,:,k) = linear_cdat(:,:,:,2*(k-1)) - ldat(:,:,:,k);
                spm_progress_bar('Set',k);
            end;
        end;
        if FirstimageType ==0;
            pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
            spm_progress_bar('Set',1);
            for k=2:perfnum,
                pdat(:,:,:,k) = cdat(:,:,:,k) - linear_ldat(:,:,:,2*(k-1));
                spm_progress_bar('Set',k);
            end;
        end;
        
    end;
    
    
    % do the sinc interpolation...
    if SubtractionType==2,
        %      sinc_ldat=zeros([VB(1).dim(1:3) 2*perfnum]);
        %      sinc_cdat=zeros([VB(1).dim(1:3) 2*perfnum]);
        sinc_ldat=zeros([VO(1).dim(1:3) 2*perfnum]);
        sinc_cdat=zeros([VO(1).dim(1:3) 2*perfnum]);
        for x=1:V(1).dim(1),
            for y=1:V(1).dim(2),
                for z=1:V(1).dim(3),
                    cdata = zeros(1,perfnum);
                    ldata = zeros(1,perfnum);
                    sinc_cdata = zeros(1,length(V));
                    sinc_ldata = zeros(1,length(V));
                    for k=1:perfnum,
                        cdata(k) = cdat(x,y,z,k);
                        ldata(k) = ldat(x,y,z,k);
                    end;
                    sincnum = fix(perfnum/Timeshift);
                    sinc_cdata=interpft(cdata,sincnum);
                    sinc_ldata=interpft(ldata,sincnum);
                    for k=1:2*perfnum,
                        sinc_cdat(x,y,z,k)= sinc_cdata(k);
                        sinc_ldat(x,y,z,k)= sinc_ldata(k);
                    end;
                end;
            end;
        end;
        
        % do the sinc subtraction....
        if FirstimageType ==0;
            for k=1:perfnum,
                pdat(:,:,:,k) = sinc_cdat(:,:,:,2*k) - ldat(:,:,:,k);
                spm_progress_bar('Set',k);
            end;
        end;
        if FirstimageType ==1;
            for k=1:perfnum,
                pdat(:,:,:,k) = cdat(:,:,:,k) - sinc_ldat(:,:,:,2*k);
                spm_progress_bar('Set',k);
            end;
        end;
    end;
    
elseif optionPCA==1
    diffs=diffs./ASLscaling;
    pdat=permute(diffs,[2 3 4 1]);
    clear diffs;
end
pdat(isnan(pdat)) = 0;
pdat(isinf(pdat)) = 0;

if PerfusionFlag == 1 && FrameFlag == 1
    % Write Perfusion image...
    for k=1:perfnum,
        VO(k) = spm_write_vol(VO(k),squeeze(pdat(:,:,:,k)));
    end;
end


% calculated the mean image...
if MeanFlag ==1 && PerfusionFlag == 1,
    Mean_dat=zeros([V(1).dim(1:3)]);
    VMP.fname = fullfile(pth,['AMean_Perf' nm(1:4) xt]);
    VMP.dt = [16 0];
    VMP.pinfo(1)=1;
    VMP = spm_create_vol(VMP);
    
    for x=1:V(1).dim(1),
        for y=1:V(1).dim(2),
            for z=1:V(1).dim(3),
                Mean_dat(x,y,z) = mean(pdat(x,y,z,:));
            end;
        end;
    end;
    
    % Write mean perfusion image...
    VMP = spm_write_vol(VMP,Mean_dat);
end;


% calculated the CBF image...
if CBFFlag ==1,
    clear ldat;
    cbfdat = zeros([VO(1).dim(1:3), perfnum]);
    cmean_dat = zeros([VO(1).dim(1:3)]);
    
    for x=1:V(1).dim(1),
        for y=1:V(1).dim(2),
            for z=1:V(1).dim(3),
                cmean_dat(x,y,z) = mean(cdat(x,y,z,:));
            end;
        end;
    end;
    clear cdat;
    
    % Write mean BOLD/Control image...
    if MeanFlag ==1 && BOLDFlag == 1,
        VMC.fname = fullfile(pth,['AMean_BOLD' nm(1:4) xt]);
        VMC.dt = [16 0];
        VMC.pinfo(1)=1;
        VMC = spm_create_vol(VMC);
        VMC = spm_write_vol(VMC,cmean_dat);
    end;
    
    for k=1:perfnum,
        for x=1:V(1).dim(1),
            for y=1:V(1).dim(2),
                for z=1:V(1).dim(3),
                    if ASLType == 3
                        Dtime = Delaytime;
                    else
                        Dtime = Delaytime + Slicetime(z);
                    end
                    if PASLModat(x,y,z) == 0
                        continue;
                    end
                    if ASLType ==0,
                        if Mask(x,y,z)<1,
                            cbfdat(x,y,z,k)=0;
                        else
                            cbfdat(x,y,z,k) = (2700*pdat(x,y,z,k)*R/alp/((exp(-Dtime*R)-exp(-(Dtime+Labeltime)*R))*PASLModat(x,y,z)));    % replacing cmean_dat with PASLModat     Danny
                        end;
                    end;
                    if ASLType ==2 || ASLType==3,
                        if Mask(x,y,z)<1,
                            cbfdat(x,y,z,k)=0;
                        else %2700 means lambda = 0.9
                            cbfdat(x,y,z,k) = (2700*pdat(x,y,z,k)*R/alp/((exp(-Dtime*R)-exp(-(Dtime+Labeltime)*R))*PASLModat(x,y,z)));    % replacing cmean_dat with PASLModat    Danny
                        end;
                    end;
                    if ASLType ==1;
                        if Mask(x,y,z)<1
                            cbfdat(x,y,z,k)=0;
                        else
                            cbfdat(x,y,z,k) = (2700*pdat(x,y,z,k)/Labeltime/alp/(exp(-(Dtime+Labeltime)*R)*PASLModat(x,y,z)));
                        end;
                    end;
                end;
            end;
        end;
        
        if FrameFlag == 1
            % Write CBF images...
            VCBF(k) = spm_write_vol(VCBF(k),cbfdat(:,:,:,k));
            spm_progress_bar('Set',k);
        end
        
    end;
    
    
    if MeanFlag ==1,
        Mean_cbfdat=zeros([VO(1).dim(1:3)]);
        Mean_cbfdat2=zeros([VO(1).dim(1:3)]);
        Mean_cbfdat3=zeros([VO(1).dim(1:3)]);
        Mean_cbfdat4=zeros([VO(1).dim(1:3)]);
        VMCBF.fname = fullfile(pth,['AMean_CBF' nm(1:4) xt]);
        VMCBF.dt = [16 0];
        VMCBF.pinfo(1)=1;%tempMo.hdr.dime.scl_slope;
        VMCBF = spm_create_vol(VMCBF);
        fileCBF=VMCBF.fname;
        
        voxelnum=0;
        vox3=0;
        vox4=0;
        zeronum=0;
        globalCBF=0;
        globalCBF2=0;
        globalCBF3=0;
        globalCBF4=0;
        globalCBFL=0;
        globalCBFR=0;
        voxL=0;
        voxR=0;
        meancontrol=0;
        halfx=round(V(1).dim(1)/2);        
        
        for x=1:V(1).dim(1),
            for y=1:V(1).dim(2),
                for z=1:V(1).dim(3),
                    Mean_cbfdat(x,y,z) = mean(cbfdat(x,y,z,:));
                    Mean_cbfdat2(x,y,z) = mean(abs(cbfdat(x,y,z,:)));
                    clear m0 n0;
                    kk0=0;
                    kk1=0;
                    for k=1:perfnum,
                        if (cbfdat(x,y,z,k) > 0.0)
                            kk1=kk1+1;
                            m0(kk1)=cbfdat(x,y,z,k);
                        else
                            kk0=kk0+1;
                            n0(kk0)=cbfdat(x,y,z,k);
                        end;
                        
                    end;
                    if kk1>0 Mean_cbfdat3(x,y,z) = mean(m0(:)); end;
                    if kk0>0 Mean_cbfdat4(x,y,z) = mean(n0(:)); end;
                    if Mean_cbfdat(x,y,z) ==0 | cmean_dat(x,y,z)< 0.1,
                        zeronum = zeronum+1;
                    else
                        voxelnum = voxelnum+1;
                        globalCBF = globalCBF+Mean_cbfdat(x,y,z);
                        globalCBF2 = globalCBF2+Mean_cbfdat2(x,y,z);
                        globalCBF3 = globalCBF3+Mean_cbfdat3(x,y,z);
                        globalCBF4 = globalCBF4+Mean_cbfdat4(x,y,z);
                        meancontrol = meancontrol+cmean_dat(x,y,z);
                        if x<halfx
                            globalCBFR = globalCBFR+Mean_cbfdat(x,y,z);
                            voxR=voxR+1;
                        else
                            globalCBFL = globalCBFL+Mean_cbfdat(x,y,z);
                            voxL=voxL+1;
                        end;
                    end;
                end;
            end;
            spm_progress_bar('Set',x*0.5);        
        end;
        
        
        mCBF=Mean_cbfdat(:); MASKm=Mask(:);
        mCBF=mean(mCBF(find(MASKm)));
        
        globalCBF = globalCBF/voxelnum;
        globalCBF2 = globalCBF2/voxelnum;
        globalCBF3 = globalCBF3/voxelnum;
        globalCBF4 = globalCBF4/voxelnum;
        meancontrol = meancontrol/voxelnum;
        
        globalCBFR = globalCBFR/voxR;
        globalCBFL = globalCBFL/voxL;
        % Write mean CBf image...
        VMCBF = spm_write_vol(VMCBF, Mean_cbfdat);
    else
        mCBF = [];
    end;
 
end;

if MeanFlag == 1
    if CBFFlag == 1
        gcbf = spm_global(VMCBF);
    end
    if BOLDFlag == 1
        gbold = spm_global(VMC);
    end
else
    gcbf = [];
    gbold = [];
end
% 
% if PerfusionFlag == 1 && FrameFlag == 1
% fprintf('\n\t Perfusion images written to: ''%s'' to \n',VO(1).fname)
% fprintf('\t  ''%s''. \n',VO(perfnum).fname)
% end
% if MeanFlag ==1 && PerfusionFlag == 1,
%     fprintf('\t Mean_perf image written to: ''%s''.\n\n',VMP.fname)
% end;
% if BOLDFlag == 1
%     if FrameFlag == 1
%         fprintf('\t BOLD images written to: ''%s''  to \n',VB(1).fname)
%         fprintf('\t ''%s'' . \n',VB(perfnum).fname)
%     end
%     if MeanFlag ==1,
%         fprintf('\t Mean_BOLD image written to: ''%s''.\n\n',VMC.fname)
%     end;
% end
% 
% if CBFFlag ==1,
%     if FrameFlag == 1
%         fprintf('\t Quantified CBF images written to: ''%s'' to \n',VCBF(1).fname)
%         fprintf('\t ''%s'' . \n',VCBF(perfnum).fname)
%     end
%     if MeanFlag == 1
%         fprintf('\t Mean Quantified CBF image written to: ''%s'' \n\n',VMCBF.fname)
%         
%         fprintf('\t the calculated voxel number is:')
%         fprintf('\t\t %8.1f\n',voxelnum)
%         fprintf('\t the zero number is:')
%         fprintf('\t\t\t %8.1f\n',zeronum)
%     end
%     if BOLDFlag == 1 && MeanFlag == 1
%         fprintf('\t the global mean BOLD control signal is:')
%         fprintf('\t\t %6.2f \n\n',meancontrol)
%     end
%     
%     fprintf('\t ---------------------------------------------------------------------\n')
%     if ASLType ==0
%         fprintf('\t -------continuous ASL [CASL]-----------------------------\n')
%     end;
%     if ASLType ==1
%         fprintf('\t -------pulsed ASL [PASL]---------------------------------\n')
%     end;
%     if ASLType ==2
%         fprintf('\t -------pseudoCASL [pCASL]--------------------------------\n')
%     end;
%     if ASLType ==3
%         fprintf('\t -------3D pseudoCASL [pCASL]-----------------------------\n')
%     end;
%     fprintf('\t ---------------------------------------------------------------------\n\n')
%     if MeanFlag == 1
%         fprintf('\t the global mean CBF signal is     :')
%         fprintf('\t %6.3f \t [ml/100g/min] \n\n',mCBF)
%         fprintf('\t ---------------------------------------------------------------------\n')
%         fprintf('\t the global    CBF>0 signal is         :')
%         fprintf('\t %6.3f',globalCBF3)
%         fprintf('\t [ml/100g/min] \n')
%         fprintf('\t the global    CBF<0 signal is         :')
%         fprintf('\t %6.3f',globalCBF4)
%         fprintf('\t [ml/100g/min] \n\n')
%     end
%     
%     
%     
%     
%     fid=fopen(prefix,'w');
%     fprintf(fid,' \n BIOCLINICA ASL BIOGEN trial \n \n');
%     fprintf(fid,' \n K. Jann @ University of Southern California, Los Angeles (version 1.1 Dez 2016) \n \n');
%     fprintf(fid,' \n File created:    %s \n \n', datestr(clock));
%     fprintf(fid,' \n TXT File stored under :    %s \n \n',prefix);
%     
%     if FrameFlag == 1
%         fprintf(fid,'\t Quantified CBF images written to: ''%s'' \n',VCBF(1).fname);
%         fprintf('\t ''%s'' . \n',VCBF(perfnum).fname)
%     end
%     if MeanFlag == 1
%         fprintf(fid,'\t Mean Quantified CBF image written to: ''%s'' \n\n',VMCBF.fname);
%         fprintf(fid,'\t the calculated voxel number is:');
%         fprintf(fid,'\t\t %8.1f\n',voxelnum);
%         fprintf(fid,'\t the zero number is:');
%         fprintf(fid,'\t\t\t %8.1f\n',zeronum);
%     end
%     if BOLDFlag == 1
%         fprintf(fid,'\t the spm global mean BOLD control signal is:');
%         fprintf(fid,'\t %6.2f \n',gbold);
%     end
%     
%     if BOLDFlag == 1 && MeanFlag == 1
%         fprintf(fid,'\t the global mean BOLD control signal is:');
%         fprintf(fid,'\t\t %6.2f \n\n',meancontrol);
%     end
%     if MeanFlag == 1
%         fprintf(fid,'\t the absolute global mean CBF signal is:');
%         fprintf(fid,'\t\t %6.3f',globalCBF2);
%         fprintf(fid,'\t [ml/100g/min] \n\n');
%     end
%     fprintf(fid,'\t ---------------------------------------------------------------------\n');
%     if ASLType ==0
%         fprintf(fid,'\t -------continuous ASL [CASL]---------------------APN (9/2009)--------\n');
%     end;
%     if ASLType ==1
%         fprintf(fid,'\t -------pulsed ASL [PASL]-------------------------APN (9/2009)--------\n');
%     end;
%     if ASLType ==2
%         fprintf(fid,'\t -------pseudoCASL [pCASL]------------------------APN (9/2009)--------\n');
%     end;
%     if ASLType ==3
%         fprintf(fid,'\t -------3D pseudoCASL [pCASL]------------------------APN (9/2009)--------\n');
%     end;
%     fprintf(fid,'\t ---------------------------------------------------------------------\n\n');
%     
%     if MeanFlag == 1
%         fprintf(fid,'\t the spm global mean CBF signal is     :');
%         fprintf(fid,'\t %6.3f \t [ml/100g/min] \n\n',gcbf);
%         fprintf(fid,'\t ---------------------------------------------------------------------\n');
%         fprintf(fid,'\t the global    CBF>0 signal is         :');
%         fprintf(fid,'\t %6.3f',globalCBF3);
%         fprintf(fid,'\t [ml/100g/min] \n');
%         fprintf(fid,'\t the global    CBF<0 signal is         :');
%         fprintf(fid,'\t %6.3f',globalCBF4);
%         fprintf(fid,'\t [ml/100g/min] \n\n');
%         fprintf(fid,'\t the global    CBF Rigth hemisphere    :');
%         fprintf(fid,'\t %6.3f',globalCBFR);
%         fprintf(fid,'\t [ml/100g/min] \n');
%         fprintf(fid,'\t the global    CBF Left  hemisphere    :');
%         fprintf(fid,'\t %6.3f',globalCBFL);
%         fprintf(fid,'\t [ml/100g/min] \n\n');
%         fprintf(fid,'\t the global mean CBF   :');
%         fprintf(fid,'\t %6.3f',mCBF);
%         fprintf(fid,'\t [ml/100g/min] \n\n');
%     end
%     
%     fclose(fid);
%     
% end;

t=toc;
%fprintf('...... CBF computing done (%5.2f sec) .\n\n',t);
spm('Pointer');
