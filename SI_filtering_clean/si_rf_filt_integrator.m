function [ unitary,uni_hist,hist_time ] = si_rf_filt_integrator(opt_params,filterType,cutoffFreqMHz)
% filterType:
%   1 - Ideal waveform
%   2 - Apply Transfer Function to the waveform
%   3 - Apply Transfer Function and then compensate over a range up to cutoffFreqMHz

%% Load relevant fields and initialize constants for Hamiltonian calculation. 

    % MY NOTES: waveform total time and timesteps
    %time_total = 4e-5; % in s
    time_total = opt_params.tot_time;
    int_timestep = 1e-8;
    % hbar = 1.054e-34;
    % hbar = 1;

    I = eye(16);

    fup = 4;
    fdown = 3;
    
    dim = 2*(fup+fdown+1);
    dim_up = 2*fup+1;
    dim_down = 2*fdown+1;
    
    % MY NOTES: create angular momentum operator matrices for the subbases
    % in the F=4 and F=3 manifolds
    upang = bgrape_make_ang_mom(fup);
    downang = bgrape_make_ang_mom(fdown);
    
    % MY NOTES: set field frequencies, phases, and amplitudes
    rf_freq = opt_params.rf_freq;
    rf_amp_x = opt_params.rf_amp_x;
    rf_amp_y = opt_params.rf_amp_y;
    mw_amp = opt_params.mw_amp;
    rf_det = opt_params.rf_det;
    
    rf_bias = rf_freq - rf_det;
    grel = -1.0032;
    hf_freq = 2*pi*(9.19263e9);
    freq_44_33_0det = hf_freq - (7*grel*rf_freq^2*(1/hf_freq)) + ((4-3*grel)*rf_freq);
    mw_det = freq_44_33_0det - ( hf_freq - (7*grel*rf_bias^2*(1/hf_freq)) + ((4-3*grel)*rf_bias) );
    mw_freq = freq_44_33_0det;

    phix = opt_params.control_fields(:,1);
    phiy = opt_params.control_fields(:,2);
    phimw = opt_params.control_fields(:,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%MAKE ALL OF THE OPERATORS IN THE HAMILTONIAN
    %%% Use the standard convention, (4,3,...,-3,-4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make the 16x16 angular momentum operators
    fx_up = [upang.jx,zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),zeros(dim_down,dim_down)];
    fy_up = -[upang.jy,zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),zeros(dim_down,dim_down)];
    fz_up = -[upang.jz,zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),zeros(dim_down,dim_down)];
    fx_down = [zeros(dim_up,dim_up),zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),downang.jx];
    fy_down = -[zeros(dim_up,dim_up),zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),downang.jy];
    fz_down = -[zeros(dim_up,dim_up),zeros(dim_up,dim_down);
        zeros(dim_down,dim_up),downang.jz];
    
    %Make the 16x16 pauli matrices for microwave transition
    mw_sigma_x = zeros(dim,dim); mw_sigma_y = zeros(dim,dim);
    mw_sigma_x(1,10) = 1; mw_sigma_x(10,1) = 1; % for |3,3> -> |4,4>
    mw_sigma_y(1,10) = i; mw_sigma_y(10,1) = -i;

    % MY NOTES: create projection matrices onto the spin manifolds
    up_proj = zeros(dim,dim); down_proj = zeros(dim,dim);
    for jj = 1:9
        up_proj(jj,jj) = 1;
    end
    for jj = 10:16
        down_proj(jj,jj) = 1;
    end
    
    
%% exact Hamiltonian     

%     H_0 = @(t) (hf_freq/2 - 16*grel*rf_bias^2*(1/hf_freq))*(up_proj - down_proj) ...
%         + rf_bias*(fz_up + grel*fz_down) + (grel*rf_bias^2 / hf_freq)*(fz_up^2 - fz_down^2);
  
%     H_RF = @(t) rf_amp_x*(cos(rf_freq*t - phix(ceil(t/opt_params.samp_time)))*(fx_up + grel*fx_down)) ...
%         + rf_amp_y*(cos(rf_freq*t - phiy(ceil(t/opt_params.samp_time)))*(fy_up + grel*fy_down));

%     H_UW = @(t) mw_amp*(cos(mw_freq*t - phimw(ceil(t/opt_params.samp_time)))) * ...
%         mw_sigma_x;

%% uW RWA Hamiltonian (BE Anderson Thesis 2013)

    H_0 =  @(t) (3*rf_bias/2 * (1+grel) - 25*grel*rf_bias^2*(1/(2*hf_freq)) ... 
        - 1/2*(mw_det - 7*rf_det))*(up_proj - down_proj) + rf_bias*(1+grel)*fz_down ...
        + (grel*rf_bias^2 / hf_freq)*(fz_up^2 - fz_down^2) - rf_det*(fz_up - fz_down);

    %H_RF in the lab frame
    
    % Prepare bandwidth-limited waveforms
    
%     omega_rf = 2*pi*1e6;
    omega_rf = opt_params.rf_freq;
%     tmax = opt_params.tot_time; %Unpadded
    tmax = opt_params.tot_time*2; %Zero-padded
    
    a = 0;
    maxdt = 1e-8;
    while tmax/(2^a) > maxdt
        a = a+1;
    end

    fprintf(strcat('The time is broken up into 2^',num2str(a),' steps\n'));

    time = linspace(0,tmax,2^a);
    dt = time(2)-time(1);
    Fs = 1/dt;

    len = length(time);
    signalX = zeros(1,len);
    signalY = zeros(1,len);

    freq = linspace(-Fs/2,Fs/2,len);
    df = freq(2)-freq(1);
    
    
    for a = 1:1:len/2
        signalX(a) = exp(1i*((omega_rf*(time(a)+0)) - opt_params.control_fields(ceil(a/((len/2)/opt_params.timesteps)),1)));
        signalY(a) = exp(1i*((omega_rf*(time(a)+0)) - opt_params.control_fields(ceil(a/((len/2)/opt_params.timesteps)),2)));
    end

    signalXInit = signalX;
    signalYInit = signalY;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%GENERATE THE TRANSFER FUNCTIONS THEN APPLY THEM TO THE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Use X_new and Y_new for latest measurements. Use X0 and Y0 for no transfer
    % function compensation.
%     [transFuncFinalX,transFuncFinalY] = makeTransferFunctions_new( './TransferFunctionX0.mat',...
%                                                                 './TransferFunctionY0.mat',...
%                                                                 freq,0,0 );
%     
%     freq_signalX = fftshift(fft(signalX));
%     freq_signalY = fftshift(fft(signalY));
% 
%     signalX = real(ifft(ifftshift(freq_signalX.*transFuncFinalX)));
%     signalY = real(ifft(ifftshift(freq_signalY.*transFuncFinalY)));
% 
%     % Find time delay
%     [convXMax,convXIndex] = max(conv(real(signalX),real(signalXInit(len:-1:1))));
%     [convYMax,convYIndex] = max(conv(real(signalY),real(signalYInit(len:-1:1))));
% 
%     timeDelayX = abs(len-convXIndex)*dt;
%     timeDelayY = abs(len-convYIndex)*dt;
% 
%     fprintf(strcat('Time Delay X-Coil: ',num2str(timeDelayX),'\n'));
%     fprintf(strcat('Time Delay Y-Coil: ',num2str(timeDelayY),'\n'));
% 
% %     for a = 1:1:len/2
% %         signalX(a) = exp(1i*((omega_rf*(time(a)+0)) - opt_params.control_fields(ceil(a/((len/2)/opt_params.timesteps)),1)));
% %         signalY(a) = exp(1i*((omega_rf*(time(a)+0)) - opt_params.control_fields(ceil(a/((len/2)/opt_params.timesteps)),2)));
% %     end
% 
%     signalX = signalXInit;
%     signalY = signalYInit;

    if filterType == 1
        timeDelayX = 0;
        timeDelayY = 0;
        [transFuncFinalX,transFuncFinalY] = makeTransferFunctions_new( './TransferFunctionX0.mat',...
                                                                './TransferFunctionY0.mat',...
                                                                freq,timeDelayX,timeDelayY);
    else
        load('./TransferFunctionX_new.mat');
        timeDelayX = TransferFunctionX(20,3)/(2*pi*1e6);
        load('./TransferFunctionY_new.mat');
        timeDelayY = TransferFunctionY(20,3)/(2*pi*1e6);
        [transFuncFinalX,transFuncFinalY] = makeTransferFunctions_new( './TransferFunctionX_new.mat',...
                                                                './TransferFunctionY_new.mat',...
                                                                freq,timeDelayX,timeDelayY);
    end

    freq_signalX = fftshift(fft(signalX));
    freq_signalY = fftshift(fft(signalY));

    signalX = real(ifft(ifftshift(freq_signalX.*transFuncFinalX)));
    signalY = real(ifft(ifftshift(freq_signalY.*transFuncFinalY)));

    % Set filter order and parameters
%     bandwidth = (1/(cutoffFreqMHz*1e6)); %Enter as 1/frequency (units of seconds)
%     nCutoffFreq = (dt*2)/bandwidth;
    
    %Butterworth Cutoff
%     butterOrder = 1;
%     [ButterA ButterB] = butter(butterOrder,nCutoffFreq);
    
    %Top-Hat Cutoff
    
    if filterType == 3
        topHat = ones(1,len);
        for a = 1:len
            if abs(freq(a)) > cutoffFreqMHz*1e6
                topHat(a) = 0;
            end
        end
        
        passBand = topHat;
        
        inverseTransFuncX = passBand./transFuncFinalX;
        inverseTransFuncY = passBand./transFuncFinalY;

        inverseTransFuncX(isinf(inverseTransFuncX)) = 1;
        inverseTransFuncY(isinf(inverseTransFuncY)) = 1;

        inverseTransFuncX(isnan(inverseTransFuncX)) = 1;
        inverseTransFuncY(isnan(inverseTransFuncY)) = 1;
        
        freq_signalX = fftshift(fft(signalX));
        freq_signalY = fftshift(fft(signalY));
        
        signalX = real(ifft(ifftshift(freq_signalX.*inverseTransFuncX)));
        signalY = real(ifft(ifftshift(freq_signalY.*inverseTransFuncY)));
    end
    
    %Xstd = std(signalX);
    %Ystd = std(signalY);
    
    %signalX = filter(ButterA,ButterB,signalX);
    %signalY = filter(ButterA,ButterB,signalY);
    
    %XModstd = std(signalX);
    %YModstd = std(signalY);

    %signalX = signalX*(Xstd/XModstd);
    %signalY = signalY*(Ystd/YModstd);

    waveform_x = @(t) interp1(time,signalX,t);
    waveform_y = @(t) interp1(time,signalY,t);
    
    H_RF_LabFrame = @(t) (rf_amp_x * waveform_x(t) * (fx_up + grel*fx_down)) + ...
                            (rf_amp_y * waveform_y(t) * (fy_up + grel*fy_down));
    
    % The transformation to the rotating frame (U_UW and U_RF)
    U_RF = @(t) expm(-1i*rf_freq*t*(fz_up - fz_down));
    U_UW = @(t) expm(-1i*(t/2)*(mw_freq-((3+4)*rf_freq))*(up_proj-down_proj));
    
    U_RWF = @(t) U_RF(t)*U_UW(t);
    
    %H_RF in the rotating frame
    H_RF = @(t) ctranspose(U_RWF(t))*H_RF_LabFrame(t)*U_RWF(t);
    
    % Microwave Hamiltonian
    H_UW = @(t) mw_amp/2*(cos(phimw(ceil(t/opt_params.samp_time))) * mw_sigma_x + ...
        sin(phimw(ceil(t/opt_params.samp_time))) * mw_sigma_y) + ...
        mw_amp/2 * si_calc_uw_RWA(t,rf_freq,phimw(ceil(t/opt_params.samp_time)));
    
    
    hammy = @(tt) ( H_0(tt) + H_RF(tt) + H_UW(tt) );

%     hammy((4e-5)-1e-15)

%% iterative brute-force quadrature numerical integration

    % MY NOTES: integrate the calculated hamiltonian to obtain the
    % effective unitary transformation through all the phase steps
    tic; 

    %Attempt at using ode45
    [unitary,uni_hist,hist_time] = unitaryEvolutionTotal_MATLAB(hammy,opt_params);
    

    
%     save('./qkt_k2_p0p99_4_40/res100_bdw1_lin_interp_rkSolved_unitaryTotal_MATLAB.mat','unitary');
%     save('./qkt_k2_p0p99_4_40/res100_bdw1_lin_interp_rkSolved_unitaryTotal_diag_6.mat','unitary_diagnostic');
    
    toc

end
