function [] = ButterworthTest(robustness,sequence,concatIndex)

% This script runs si_calc_filt_opt_integrator_proper for a given
% opt_params and finds when bandwidth limitations start to
% break down the fidelity of the resulting unitary.

%QKT waveforms
% load('./qkt_k2_p0p99_4_40/varyingKappa2kap_4step_40tot_n1_s79.mat');
% load('./refinedFields3/kappa2_samp4_tot40_jobnum139_3_new.mat');

%Benchmark waveforms
% load('../waveforms/concat_unibench_20130715_R_d16/concat_unibench_20130715_R_1_3.mat');
load('../waveforms/concat_unibench_20130715_NR_d16/concat_unibench_20130715_NR_1_3.mat');

%Load Time delays
load('./ButterworthFilterTimeDelays.mat');

if robustness == 0
   robString = 'NR';
else
   robString = 'R';
end

load(strcat('../waveforms/concat_unibench_20130715_',robString,...
'_d16/concat_unibench_20130715_',robString,'_',num2str(sequence),...
'_',num2str(concatIndex),'.mat'));

res = 1;

%timeDelay = linspace(2e-7,4e-7,res);
cutoffFreqMHz = linspace(25,25,res);

outFid = zeros(res,1);
uniEigFid = zeros(res,1);
uniFid = zeros(res,1);

unitary = zeros(16,16,res);

for a = 1:res
    
    [unitary(:,:,a),uni_hist,hist_time] = si_calc_filt_opt_integrator_proper(opt_params,TimeDelay(25),cutoffFreqMHz(a));
    
    [outFid(a),uniEigFid(a),uniFid(a)] = checkUniFid_fn(unitary(:,:,a),opt_params.uni_final);
  
    save(strcat('./BandwidthRuns/ButterworthRuns/uniHist_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'_',num2str(cutoffFreqMHz(a)),'MHz.mat'),'uni_hist');
    save(strcat('./BandwidthRuns/ButterworthRuns/histTime_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'_',num2str(cutoffFreqMHz(a)),'MHz.mat'),'hist_time');

end

%unitary = si_calc_filt_opt_integrator_proper(opt_params,0,1);
%
%[stateFid,uniFid,cost] = checkUniFid_fn(unitary,opt_params.uni_final);


%figure(4)
%plot(timeDelay,uniFid,'LineWidth',2);
%xlabel('Time Delay (s)');
%ylabel('Unitary Fidelity');
%
%figure(5)
%plot(timeDelay,stateFid,'LineWidth',2);
%xlabel('Time Delay (s)');
%ylabel('State Fidelity');
%
%figure(6)
%plot(timeDelay,cost,'LineWidth',2);
%xlabel('TimeDelay (s)');
%ylabel('Cost Function');
    
save(strcat('./BandwidthRuns/ButterworthRuns/uniFid_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'uniFid');
save(strcat('./BandwidthRuns/ButterworthRuns/unitary_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'unitary');

%save('./ToleranceTest/uniFid.mat','uniFid');
%save('./ToleranceTest/unitary.mat','unitary');

end
