function [ transFuncFinalX,transFuncFinalY ] = makeTransferFunctions_new( TransFuncX,TransFuncY,freq,timeDelayX,timeDelayY )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Load and Prepare Measured Transfer Function Parameters NEW (Picoscope Measurements)
    load(TransFuncX);
    load(TransFuncY);
    
    dataLength = length(TransferFunctionX(:,1));

    frequency = TransferFunctionX(:,1);
    xCoil_amp = TransferFunctionX(:,2);
    yCoil_amp = TransferFunctionY(:,2);

    phaseDelayX = TransferFunctionX(:,3);
    phaseDelayY = TransferFunctionY(:,3);
    
%     timeDelay = 0.25e-6;
%     timeDelay = 0;
    phaseDelayX = phaseDelayX - (timeDelayX*frequency*2000*pi);
    phaseDelayY = phaseDelayY - (timeDelayY*frequency*2000*pi);

    %See effects of inaccuracies in the Transfer Function Measurement
    randRange_Amp = 0; %0.5;
    xCoil_amp = xCoil_amp+(randRange_Amp*xCoil_amp.*rand(size(xCoil_amp)))-(xCoil_amp*randRange_Amp/2);
    yCoil_amp = yCoil_amp+(randRange_Amp*yCoil_amp.*rand(size(yCoil_amp)))-(yCoil_amp*randRange_Amp/2);

    randRange_Phase = 0; %2*pi/20;
    phaseDelayX = phaseDelayX+(randRange_Phase*rand(size(phaseDelayX)))-(randRange_Phase/2);
    phaseDelayY = phaseDelayY+(randRange_Phase*rand(size(phaseDelayY)))-(randRange_Phase/2);
    
    %% Prepare Transfer Functions

    freq_lim = length(phaseDelayX);
    frequency_kHz = frequency(1:1:freq_lim,1);

    % Set the amplitude of the spectrum at 1MHz to be equal to 1
    xCoil_amp = xCoil_amp/xCoil_amp(20);
    yCoil_amp = yCoil_amp/yCoil_amp(20);

    % transferFuncX = zeros(floor(length(freq)/2)-1,1);
    transferFuncX_pos = xCoil_amp(1:1:dataLength) .* exp(-1i*phaseDelayX); %For positive part of transfer function
    transferFuncX_neg = xCoil_amp(1:1:dataLength) .* exp(1i*phaseDelayX); %For negative part of transfer function
    % transferFuncY = zeros(floor(length(freq)/2)-1,1);
    transferFuncY_pos = yCoil_amp(1:1:dataLength) .* exp(-1i*phaseDelayY); %For positive part of transfer function
    transferFuncY_neg = yCoil_amp(1:1:dataLength) .* exp(1i*phaseDelayY); %For negative part of transfer function


    % Make the transfer function symmetric in the negative frequency regime

    symm_TF_X = zeros(1,(2*length(transferFuncX_pos))+1);
    symm_TF_X(ceil(length(symm_TF_X)/2)+1:1:length(symm_TF_X)) = transferFuncX_pos; %set positive half
    symm_TF_X(1:1:floor(length(symm_TF_X)/2)) = transferFuncX_neg(length(transferFuncX_neg):-1:1); %set negative half
    symm_TF_X(ceil(length(symm_TF_X)/2)) = 1; %set the 0 frequency value to unity

    symm_TF_Y = zeros(1,(2*length(transferFuncY_pos))+1);
    symm_TF_Y(ceil(length(symm_TF_Y)/2)+1:1:length(symm_TF_Y)) = transferFuncY_pos; %set positive half
    symm_TF_Y(1:1:floor(length(symm_TF_Y)/2)) = transferFuncY_neg(length(transferFuncY_neg):-1:1); %set negative half
    symm_TF_Y(ceil(length(symm_TF_Y)/2)) = 1; %set the 0 frequency value to unity

    % Interpolate the transfer function for finer frequency resolution

    transFuncX_interp = @(f) interp1([1:1:(2*dataLength)+1],symm_TF_X,(f/50000)+(dataLength+1));
    transFuncY_interp = @(f) interp1([1:1:(2*dataLength)+1],symm_TF_Y,(f/50000)+(dataLength+1));

    transFuncFinalX = transFuncX_interp(freq);
    transFuncFinalX(isnan(transFuncFinalX)) = 1;

    transFuncFinalY = transFuncY_interp(freq);
    transFuncFinalY(isnan(transFuncFinalY)) = 1;

end

