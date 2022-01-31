

function [globalCBF, DVARS, FD, tSNR, tSD, CBF, adjCBF] = LOFT_CBFquantify(images, Mzero, params)
%%[globalCBF, DVARS, FD, tSNR, tSD, CBF, adjCBF] = LOFT_CBFquantify(images, Mzero, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, PLD, LabelTime, Slicetime, threshold, optionPCA, ASLscaling, M0scaling)
%
%Quality Control (QC) parameter calculation for Arterial Spin Labeling MRI data
%Code written by Kay Jann, Ph.D, University of Southern California
%Code parts provided by Robert Smith, Edgar Rios Piedra, Xingfeng Xiao, Danny JJ Wang
%July 2016 (version 3.0 updated September 2017)
%F
%
%
%INPUT PARAMETERS
%
%images, image files of alternating ASL label and control images 
%MO, Zero Magnetization Images, required for all sequences. If not acquired during scan use first control image (no background suppression for M0 tho) 
%FieldStrength, scanner field strength [0=1.5T, 1=3T]
%ASLType, ASL sequence type (0=CASL, 1=PASL, 2=pCASL, 3=3D pCASL with BS)
%FirstimageType, integer variable indicating the type of first image - 0:control; 1:labeled
%SubtractionOrder, depending on first image time select Even-Odd(Img2-Img1)-> 0 |  or Odd-Even(Img1-Img2) -> 1 , so that subtraction is Control-Label 
%SubtractionType, 0 simple subtraction / 1 surround subtraction / 2 sinc-subtraction
%PLD, post-label delay time in seconds (usual ranges form 1sec to 2sec for single delay ASL sequences
%LabelTime, duration of labeling pulse/pulses in seconds
%Slicetime, slice-timing in 2D acquisitions, i.e. time between two slice readouts [in milliseconds] / set to 0 in case of single shot aquisitions
%threshold, threshold for mask (percentage of maximum intensity of M0 image, suggested value is 0.1).
%optionPCA, [0 1], option to perform PCA decomposition on ASL images before estimating CBF
%ASLscaling, MOscaling, usually 1 for both, except for Philips sequences, scaling factor can be found in dicom header as Private_2005_100e for ASL and M0 scans respectively

FieldStrength = params.FieldStrength;
ASLType = params.ASLType;
FirstimageType = params.FirstimageType;
SubtractionOrder = params.SubtractionOrder;
SubtractionType = params.SubtractionType;
PLD = params.PLD;
LabelTime = params.LabelTime;
Slicetime = params.Slicetime;
threshold = params.threshold;
optionPCA = params.optionPCA;
ASLscaling = params.ASLscaling;
M0scaling = params.M0scaling;

ThreshFlag = params.ThresholdFlag;
FrameFlag = params.FrameFlag;
PerfusionFlag = params.PerfusionFlag;
CBFFlag = params.CBFFlag;
MeanFlag = params.MeanFlag;
BOLDFlag = params.BOLDFlag;
AnatFlag = params.AnatFlag;
BSScaling = params.BSScaling;
%% Main Program

number_of_files=numel(images);
for i=1:number_of_files
    [PATHSTR,NAME,EXT] = fileparts(images{i});
    imagesr(i,:)=strcat(PATHSTR, '/r', NAME, EXT);
    imagesraw(i,:)=strcat(PATHSTR, '/', NAME, EXT);
end

% separate label and control images for motion correction 
if FirstimageType==0 %first image is control
    u=1; 
    for i=1:2:number_of_files; 
        [PATHSTR,NAME,EXT] = fileparts(images{i});
        path_ctr(u,:) = PATHSTR;
        name_ctr(u,:) = strcat(PATHSTR, '/', NAME,  EXT, ',1');
        u=u+1; 
    end
    u=1;  
    for i=2:2:number_of_files; 
        [PATHSTR,NAME,EXT] = fileparts(images{i});
        path_lbl(u,:) = PATHSTR;
        name_lbl(u,:) = strcat(PATHSTR, '/', NAME,  EXT, ',1');   % NEW - same for label
        u=u+1; 
    end
    name_lbl=cellstr(name_lbl(:,:));
elseif FirstimageType==1 %first image is label
    u=1; 
    for i=1:2:number_of_files; 
        [PATHSTR,NAME,EXT] = fileparts(images{i});
        path_lbl(u,:) = PATHSTR;
        name_lbl(u,:) = {[PATHSTR, '/', NAME,  EXT, ',1']};   % NEW - will contain the complete address of the control images once they are resliced, cool for input in 3dgrase
        u=u+1; 
    end
    u=1;  
    for i=2:2:number_of_files; 
        [PATHSTR,NAME,EXT] = fileparts(images{i});
        path_ctr(u,:) = PATHSTR;
        name_ctr(u,:) = {[PATHSTR , '/', NAME,  EXT, ',1']};   % NEW - same for label
        u=u+1; 
    end
end

name_ctr=cellstr(name_ctr(:,:));
name_lbl=cellstr(name_lbl(:,:));

ctrx = name_ctr(:,:);
ctrx = cellstr(ctrx);
lblx = name_lbl(:,:);
lblx = cellstr(lblx);

if iscell(Mzero)
    [PATHSTR,NAME,EXT] = fileparts(Mzero{1});
    clear Mzero
        path_Mzero(1,:) = PATHSTR;
        Mzero(1,:) = strcat(PATHSTR, '/', NAME,  EXT);
        %display(Mzero)
        
end

temp=load_untouch_nii(imagesr(1,:)); 
 
V=spm_vol(deblank(imagesr)); %changed to imagesr from images

matrix = zeros([V(1).dim,number_of_files]);
 for i=1:number_of_files
    matrix(:,:,:,i) = spm_read_vols(V(i))./temp.hdr.dime.scl_slope;
 end
 matrix(isnan(matrix)) = 0;
 matrix(isinf(matrix)) = 0;

% PCA COMPUTATION
if optionPCA
    diffs = PCA(matrix,FirstimageType,SubtractionOrder); %output in 4D T X Y Z
end

% CBF ESTIMATION
rawASL=char(cellstr(imagesr));
if ASLType==0 || ASLType==1 || ASLType==2 || ASLType==3
    %display('running 2D')
    if optionPCA==0
        [gcbf, cbfdat, mmCBF, Mask] = perf_reconstruct_single_delay(rawASL, FieldStrength, ASLType, BSScaling,FirstimageType, SubtractionOrder, SubtractionType, ThreshFlag, threshold, FrameFlag, PerfusionFlag, CBFFlag, MeanFlag, BOLDFlag,AnatFlag, PLD, LabelTime, Slicetime, Mzero, optionPCA, ASLscaling, M0scaling);
    elseif optionPCA==1
        [gcbf, cbfdat, mmCBF, Mask] = perf_reconstruct_single_delay(rawASL, FieldStrength, ASLType, BSScaling,FirstimageType, SubtractionOrder, SubtractionType, ThreshFlag, threshold, FrameFlag, PerfusionFlag, CBFFlag, MeanFlag, BOLDFlag,AnatFlag, PLD, LabelTime, Slicetime, Mzero, optionPCA, ASLscaling, M0scaling, diffs);
    end
    
elseif ASLType==4
    %display('running 3D') 
%     if optionPCA==0
%         [gcbf, cbfdat, mmCBF, Mask] = perf_reconstruct_3DGRASE(rawASL, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, 1, threshold, CBFFlag, MeanFlag, PLD, LabelTime, Slicetime, 1650, Mzero, optionPCA, ASLscaling, M0scaling);
%     elseif optionPCA==1
%         [gcbf, cbfdat, mmCBF, Mask] = perf_reconstruct_3DGRASE(rawASL, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, 1, threshold, CBFFlag, MeanFlag, PLD, LabelTime, Slicetime, 1650, Mzero, optionPCA, ASLscaling, M0scaling, diffs);
%     end
    error('This toolbox does not support multi-delay ASL. Please use multi-delay ASL toolbox instead!');
end


% DVARS CALCULATION 
try
    load('./perf/meanFD.mat');
    if FirstimageType==1 %first image is control
        DVARS.control = GetDVARS(images(1:2:number_of_files),number_of_files/2,Mask); 
        DVARS.label = GetDVARS(images(2:2:number_of_files),number_of_files/2,Mask);
        FD.label=meanFD.odd;FD.control=meanFD.even;
    elseif FirstimageType==0 %first image is control
        DVARS.label = GetDVARS(images(1:2:number_of_files),number_of_files/2,Mask); 
        DVARS.control = GetDVARS(images(2:2:number_of_files),number_of_files/2,Mask);
        FD.label=meanFD.even;FD.control=meanFD.odd;
    end
    DVARS.label=mean(DVARS.label,2);DVARS.control=mean(DVARS.control,2);
catch
    warning('No meanFD file was found in perf folder!');
    FD.label = 0;
    FD.control = 0;
end


% tSNR CALCULATION
if FirstimageType==1 %first image is control
    [tSNR.control,tSD.control] = calculateSNR(matrix(:,:,:,1:2:number_of_files),Mask); 
    [tSNR.label,tSD.label] = calculateSNR(matrix(:,:,:,2:2:number_of_files),Mask); 
elseif FirstimageType==0 %first image is control
    [tSNR.label,tSD.label] = calculateSNR(matrix(:,:,:,1:2:number_of_files),Mask); 
    [tSNR.control,tSD.control] = calculateSNR(matrix(:,:,:,2:2:number_of_files),Mask); 
end
tSNR.control=mean(tSNR.control);tSNR.label=mean(tSNR.label);
tSD.control=mean(tSD.control);tSD.label=mean(tSD.label);

[tSNR.CBF, tSD.CBF]=calculateSNR(cbfdat,Mask);
tSNR.CBF=mean(tSNR.CBF);
tSD.CBF=mean(tSD.CBF);


globalCBF=gcbf;
mCBF=mean(cbfdat,4); lmCBF=mCBF(:); CBF=mean(lmCBF(find(Mask(:))));
lowCBF=sum(find(lmCBF(logical(Mask(:))<-20)));
highCBF=sum(find(lmCBF(logical(Mask(:))>200)));
adjCBF=lmCBF(logical(Mask(:)));
adjCBF=adjCBF(find(adjCBF>20)); adjCBF=adjCBF(find(adjCBF<150)); 
adjCBF=mean(adjCBF);
CBF=CBF;

% path1=pwd;
% cd(pwd);
% [q prefix]=fileparts(images{1});
% fid=fopen(['QC_' prefix '.txt' ],'w');
%      fprintf(fid,'\n LOFT lab CBF toolbox \n \n');
%      fprintf(fid,'\n K. Jann @ University of Southern California, Los Angeles (version 1.1 Dez 2016) \n \n');
%      fprintf(fid,'\n File created:    %s \n \n', datestr(clock));
%      fprintf(fid,'\n TXT File stored under :    %s \n \n', ['QC_' prefix ]);
%      fprintf(fid,'\t the DVARS for label images is:');
%      fprintf(fid,'\t %6.3f',DVARS.label);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the DVARS for control images is:');
%      fprintf(fid,'\t %6.3f',DVARS.control);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the SD for label images is:');
%      fprintf(fid,'\t\t %6.3f',tSD.label);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the SD for control images is:');
%      fprintf(fid,'\t\t %6.3f',tSD.control);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the tSNR for label images is:');
%      fprintf(fid,'\t\t %6.3f',tSNR.label);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the tSNR for control images is:');
%      fprintf(fid,'\t %6.3f',tSNR.control);
%      fprintf(fid,'\t\n\n');  
%      fprintf(fid,'\t the FD for label images is:');
%      fprintf(fid,'\t\t %6.3f',FD.label);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the FD for control images is:');
%      fprintf(fid,'\t\t %6.3f',FD.control);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the globalCBF  is:');
%      fprintf(fid,'\t\t\t %6.3f',mmCBF);
%      fprintf(fid,'\t [ml/100g/min] \n\n');
%      fprintf(fid,'\t the adjusted is:');
%      fprintf(fid,'\t\t\t %6.3f',adjCBF);
%      fprintf(fid,'\t [ml/100g/min] \n\n');
%      fprintf(fid,'\t the SD for CBF images is:');
%      fprintf(fid,'\t\t %6.3f',tSD.CBF);
%      fprintf(fid,'\t\n\n');
%      fprintf(fid,'\t the tSNR for CBF images is:');
%      fprintf(fid,'\t\t %6.3f',tSNR.CBF);
%      fprintf(fid,'\t\n\n');
%      fclose(fid);

end
%% X 
%%%%%%%%%%%%%%%%% ADDITIONAL SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%
%% FUNCTION for computing DVARS
function DVARS = GetDVARS(Filenames,NumberOfFiles,Mask)
%computes signal variation across volumes
 
   % Mask=logical(Mask);
    %Take first image to get stuff
    FileTest =  char(Filenames(1,1))  ;                     
    Image1 = load_untouch_nii(FileTest);                            
    %view_nii(Image1);
    BruinBrainTest = getfield(Image1,'img');                
    SizeBruinTest = size(BruinBrainTest);                   
    X = SizeBruinTest(1);                                
    Y = SizeBruinTest(2);                           
    Z = SizeBruinTest(3);       

    DVARS(NumberOfFiles) = 0;                           
    DVARS(1) = 0;

    %h = waitbar(0,'Please wait...');
    for w=2:NumberOfFiles

        GrabFile1 =  char(Filenames(1,w-1));              
        Image1 = load_untouch_nii(GrabFile1); 
        BruinBrain1 = getfield(Image1,'img');
        BruinBrain1=  uint16(BruinBrain1).*uint16(Mask);

        GrabFile2 =  char(Filenames(1,w)); 
        Image2 = load_untouch_nii(GrabFile2);
        BruinBrain2 = getfield(Image2,'img');
        BruinBrain2=  uint16(BruinBrain2).*uint16(Mask);

        % Here I have both brains matrix stored in BruinBrain1 and BruinBrain2

        SliceA(X,Y)=0;  ValueA=0;  ArrayA=[];  MeansA(Z)=0; 
        SliceB(X,Y)=0;  ValueB=0;  ArrayB=[];  MeansB(Z)=0;

        for i=1:Z
            
           SliceA = BruinBrain1(:,:,i); 
           SliceB = BruinBrain2(:,:,i);
           
           for x=1:X                                       
               for y=1:Y
                   ValueA = SliceA(x, y);
                   ValueB = SliceB(x, y);
                   if(ValueA>0 && ValueB>0)                 
                       ArrayA = [ArrayA(1:end) ValueA];   
                       ArrayB = [ArrayB(1:end) ValueB];    
                   end
               end
           end

           SiseA = size(ArrayA);                                    
           SiseB = size(ArrayB);                                
   
           if (SiseA(2)==0) ArrayA = [ArrayA(1:end) 0]; end     % In case the whole slice is zero,
           if (SiseB(2)==0) ArrayB = [ArrayB(1:end) 0]; end     % adding a zero to avoid NaN
           
           MeansA(i) = mean(ArrayA);
           MeansB(i) = mean(ArrayB);

           ArrayA = [];                                  
           ArrayB = [];

        end

        SliceMeanA = 0;
        SliceMeanB = 0;
        Squares(Z) = 0;                                   

        for i=1:Z

            SliceMeanA = MeansA(i);
            SliceMeanB = MeansB(i);

            Resta = (SliceMeanB - SliceMeanA);
            Square = Resta*Resta;

            Squares(i) = Square;

        end

        % Here I got an array called squares that contains the substractes slices al cuadrado - useful to get the RMS

        MeanSquares = mean(Squares);                          % Spatial average over the whole brain, requiered by formula                         

        RMS = sqrt(MeanSquares);                              % Root Mean Square
        
        DVARS(w) = RMS;
        
        % RMS is the value of the rate of change in the bold signal across the brain in this two brains
        % And DVARS just stores that value for every iteration
        
    %disp('Finished one pair of BruinBrains')
    %waitbar(w / NumberOfFiles)
    end
    %close(h)
    
    %Results
    %DVARS
    %plot(DVARS)
end

%% FUNCTION to Apply PCA to vNav_pCASL data
function diffsVol = PCA(mo_moco_orig,FirstimageType,SubtractionOrder)
%
% Input: mo_moco, interleved control/label images, size: [Nx, Ny, Nt]
% Output: pca corrected perfusion images, size: [Nt/2, Nx, Ny]
% Reshape input data to a Nx*Ny-by-Nt matric and apply PCA. (source data in mosaic layout)
% Nt is total number of time point, Nt/2 is number of "variables"


% [Nx,Ny,Nz,Nt] = size(mo_moco);
[Nx,Ny,Nz,Nt] = size(mo_moco_orig);

%for z=1:Nz 
mo_moco=squeeze(mo_moco_orig(:,:,:,:));    
% mo_moco_pre = reshape(mo_moco,Nx*Ny*Nz,Nt);
mo_moco_pre = reshape(mo_moco,Nx*Ny*Nz,Nt);
xd = size(mo_moco_pre,1);
zd = size(mo_moco_pre,2)/2;  % number of measurements


if FirstimageType==0 %first image is control
    refs  = mo_moco_pre(:,1:2:end);
    refs  = double(refs');
    if SubtractionOrder==0 %Even-Odd(Img2-Img1)
        diffs = mo_moco_pre(:,2:2:end) - mo_moco_pre(:,1:2:end);
    elseif SubtractionOrder==1 %Odd-Even(Img1-Img2)
        diffs = mo_moco_pre(:,1:2:end) - mo_moco_pre(:,2:2:end);
    end
    diffs = double(diffs');
elseif FirstimageType==1 %first image is label
    refs  = mo_moco_pre(:,2:2:end);
    refs  = double(refs');
    if SubtractionOrder==0 %Even-Odd(Img2-Img1)
        diffs = mo_moco_pre(:,2:2:end) - mo_moco_pre(:,1:2:end);
    elseif SubtractionOrder==1 %Odd-Even(Img1-Img2)
        diffs = mo_moco_pre(:,1:2:end) - mo_moco_pre(:,2:2:end);
    end
    diffs = double(diffs');
end  

% Calculate eigenvalues and eigenvectors
[q d] = eig(refs*refs');
d= diag(d);
eimages = q'*double(refs); % project control images to eigen vector to generate eigen images

% Normalize eigen images
d = sqrt(d);

for i = 1:zd
    eimages(i,:) = eimages(i,:)./d(i);
end


% Assume first 1:zd-1 vectors are "noise bases"
% Project diffusion images to "noise bases" to generate noise components of
% perfusion images.
Noise_diffs = ((diffs)*eimages(1:zd-1,:)')*(eimages(1:zd-1,:));
diffs = (diffs) - Noise_diffs;
diffs = reshape(diffs,zd,Nx,Ny,Nz);

diffsVol(:,:,:,:)=diffs;
%end

% diff  = reshape(median(diffs),Nx,Ny);
% diff_smooth = diff;
% diff_smooth(diff_smooth <= 0.05*max(max(abs(diff_smooth)))) = 0;
% GaussianKernel = fspecial('Gaussian',[5 5],0.2);
% diff_smooth = imfilter(diff_smooth,GaussianKernel);

% figure,
% subplot(1,2,1);imagesc(abs(diff),[0 5]);colormap gray;axis off;
% subplot(1,2,2);imagesc(abs(diff_smooth),[0 5]);colormap gray;axis off;

end

%% FUNTION for computing tSNR and tSD
function [tSNR, SD] = calculateSNR(matrix,Mask)
idx=find(Mask);
dim=size(matrix);
rmat=reshape(matrix,[dim(1)*dim(2)*dim(3),dim(4)]);
SD=std(double(rmat(idx,:)),0,2); %SD of timecourse
mTC=mean(rmat(idx,:),2); %mean ot timecourse
%idxnonzeromean=find(mTC~=0);
idxnonzeromean=find(SD~=0);
tSNR=mTC(idxnonzeromean)./SD(idxnonzeromean);

%tSD=mean(SD);
%tSNR=mean(tSNR);

end

%% END