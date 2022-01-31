%% FUNCTION for computing FD
function meanFD = JointMotionParameters(motionFile)
   
    %[File, pathos] = uigetfile({'*.txt'},'Select a motion parameters file');
    
    %Original_Directory = pwd;
    %cd(motionPath);
    
    TXT  = load(char(motionFile));
    %cd(Original_Directory);

    S = size(TXT);                                                        
    Size = S(1);                                                           

    %%% Building new matrix with converted rotational parameters %%%
    Parameters = zeros(size(TXT));
    Parameters(:,1) = TXT (:,1);                                           % Copy
    Parameters(:,2) = TXT (:,2);
    Parameters(:,3) = TXT (:,3);
    JoinedMotionParameters = zeros(Size,1);                               
    Alpha = 0; Beta = 0; Gamma = 0;
    NumberOfTemporalSeries = Size;

    for i=1:NumberOfTemporalSeries

        Alpha = TXT(i,4);
        Beta  = TXT(i,5);
        Gamma = TXT(i,6);
        
        % Assuming sphere of 50mm radius
        % Approx mid distance from cerebral cortex to the center of the head 

        Alpha = (2)*(pi)*(50)*(Alpha/360);
        Beta = (2)*(pi)*(50)*(Beta/360);
        Gamma = (2)*(pi)*(50)*(Gamma/360);
        
        
        Parameters(i,4) = Alpha;
        Parameters(i,5) = Beta;
        Parameters(i,6) = Gamma;

    end
   
    % Now to compute the difference between the numbers
    
    JoinedMotionParameters(1) = 0;                                         % First
    TempValue = [0 0 0 0 0 0];
    
    for i=2:Size                                                           % From the 2nd one
        for j=1:6
            Value1 = Parameters(i-1,j);
            Value2 = Parameters(i,j);
            Delta = abs(Value1-Value2);
            TempValue(j)= Delta*Delta;   % ^2
        end
        
        Suma = TempValue(1)+TempValue(2)+TempValue(3)+TempValue(4)+TempValue(5)+TempValue(6);   % Already squared
        Suma = Suma/6;
        FramewiseDisplacement = sqrt(Suma);                                % RMS value

        JoinedMotionParameters(i) = FramewiseDisplacement;     
    end
    
    meanFD = mean(JoinedMotionParameters);
    
    %plot(JoinedMotionParameters)
    
end 