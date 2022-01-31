json_file = dir('*.json');
json_file = regexpi({json_file.name},strcat(subject_ID,'_asl\.json$'),'match');
json_file = [json_file{:}];
if isempty(json_file)
    error('json file for perf data is not found (Should be subject_ID_asl.json in perf folder)!\n');
end
fid = fopen(json_file{1}); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
[ASL_info,~] = parse_json(str);
ASL_info = ASL_info{1};
params = struct();

try
    if ASL_info.MagneticFieldStrength == 1.5
        params.FieldStrength = 0;
    elseif ASL_info.MagneticFieldStrength == 3
        params.FieldStrength = 1;
    elseif ASL_info.MagneticFieldStrength == 7
        params.FieldStrength = 2;
    else
        error('Not a valid field strength!\n');
    end
catch
    params.FieldStrength = input('FieldStrength of scanner, 0:1.5T; 1:3T;3:7T?\n');
end

try
    switch ASL_info.ArterialSpinLabelingType
        case 'CASL'
            params.ASLType = 0;
        case 'PASL'
            params.ASLType = 1;
        case 'PCASL'
            if strcmp(ASL_info.MRAcquisitionType,'2D')
                params.ASLType = 2;
            elseif strcmp(ASL_info.MRAcquisitionType,'3D')
                if length(ASL_info.PostLabelingDelay) == 1
                   params.ASLType = 3;
                else
                   params.ASLType = 4;
                end
            end 
        otherwise
            error('Error, no such ASL type! Try again!\n');
    end
catch
    params.ASLType = input('ASL type, 0:CASL; 1:PASL;2:PCASL(2D);3:PCASL(3D);4:Multi-delay PCASL(3D)?\n');
end

try
    if strcmp(config_info.first_image_type,'control')
        params.FirstimageType = 0;
    elseif strcmp(config_info.first_image_type,'label')
        params.FirstimageType = 1;
    else
        error('Error, no such first image type! Try again!\n');
    end
catch
    params.FirstimageType = input('First image type, 0:control; 1:label;\n');
end

if params.FirstimageType == 0
    params.SubtractionOrder = 1;
else 
    params.SubtractionOrder = 0;
end

params.SubtractionType = 3;
if params.SubtractionType == 3
    params.optionPCA = 1;
else
    params.optionPCA = 0;
end

try
    if params.ASLType == 4
        params.PLD = cell2mat(ASL_info.PostLabelingDelay);
    else
        params.PLD = ASL_info.PostLabelingDelay;
    end
    if isempty(ASL_info.PostLabelingDelay)
        error('PLD is not provided in perf json!\n');
    end    
catch
    if params.ASLType == 4
        n_delays = input('How many PLDs in your multi-delay 3D pCASL data?;\n');
        params.PLD = input('PLD in sec ([1.5,1.75,2,...]);\n');
        if length(params.PLD)~=n_delays
            error('PLDs you entered does not match your data!');
        end
    else
        params.PLD = input('PLD in sec;\n');
    end
end

try
    params.LabelTime = ASL_info.LabelingDuration;
    if isempty(ASL_info.LabelingDuration)
        error('PLD is not provided in perf json!\n');
    end   
catch
    params.LabelTime = input('Label duration in sec;\n');
end

try
    params.Slicetime = cell2mat(ASL_info.SliceTiming);
catch
    if params.ASLType == 3 || params.ASLType == 4
        params.Slicetime = 0;
    else
        n_slices = input('How many slices in your 2D data?;\n');
        params.Slicetime = input('Slice time in sec ([1.5,1.75,2,...]);\n');
        if length(params.Slicetime)~=n_slices
            error('Slice time you entered does not match your data!\n');
        end
    end
end
params.threshold = 0.1;
params.ASLscaling = 1;
params.M0scaling = 1;

try
    params.AnatFlag = config_info.T1w_flag;
    if isempty(config_info.T1w_flag)
        error('T1w_flag is not defined in config file!');
    end
catch
    params.AnatFlag = input('Do you have T1w data, 0:no, 1: yes?;\n');
end

if params.AnatFlag == 1
    params.ThresholdFlag = 0;
else
    params.ThresholdFlag = 1;
end

try
    params.FrameFlag = config_info.frame_flag;
    if isempty(config_info.frame_flag)
        error('frame_flag is not defined in config file!');
    end
catch
    params.FrameFlag = input('Do you want to save data of each frame, 0:no, 1: yes?;\n');
end

try
    params.PerfusionFlag = config_info.perfusion_flag;
    if isempty(config_info.perfusion_flag)
        error('perfusion_flag is not defined in config file!');
    end
catch
    params.PerfusionFlag = input('Do you want to save perfusion data, 0:no, 1: yes?;\n');
end

try
    params.CBFFlag = config_info.CBF_flag;
    if isempty(config_info.CBF_flag)
        error('CBF_flag is not defined in config file!');
    end
catch
    params.CBFFlag = input('Do you want to save CBF data, 0:no, 1: yes?;\n');
end

try
    params.MeanFlag = config_info.mean_flag;
    if isempty(config_info.mean_flag)
        error('mean_flag is not defined in config file!');
    end
catch
    params.MeanFlag = input('Do you want to save Mean data, 0:no, 1: yes?;\n');
end

try
    params.BOLDFlag = config_info.BOLD_flag;
    if isempty(config_info.BOLD_flag)
        error('BOLD_flag is not defined in config file!');
    end
catch
    params.BOLDFlag = input('Do you want to save BOLD data, 0:no, 1: yes?;\n');
end

try
    if ASL_info.BackgroundSuppression
        params.BSScaling = 0.95^ASL_info.BackgroundSuppressionNumberPulses;
    else
        params.BSScaling = 1;
    end
catch
    BackgroundSuppressionNumberPulses = input('How many BS pulses were used?(0:no BS, 1: 1 BS pulse, 2: 2 BS pulses, ...)\n');
    params.BSScaling = 0.95^BackgroundSuppressionNumberPulses;
end