%function [] = SI_TransFuncTest()
function [] = SI_TransFuncTest(robustness,sequence,concatIndex,cutoff)

%cutoff must be an integer 1-22

% This script runs si_calc_filt_opt_integrator_proper for a given
% opt_params

%QKT waveforms
% load('./qkt_k2_p0p99_4_40/varyingKappa2kap_4step_40tot_n1_s79.mat');
% load('./refinedFields3/kappa2_samp4_tot40_jobnum139_3_new.mat');
% load('./varyingKappa2kap_4step_40tot_n1_s79.mat');

%Benchmark waveforms
% load('../waveforms/concat_unibench_20130715_R_d16/concat_unibench_20130715_R_1_3.mat');
% load('../waveforms/concat_unibench_20130715_NR_d16/concat_unibench_20130715_NR_1_3.mat');

if robustness == 0
   robString = 'NR';
else
   robString = 'R';
end

%cluster load path
% load(strcat('../waveforms/concat_unibench_20130715_',robString,...
% '_d16/concat_unibench_20130715_',robString,'_',num2str(sequence),...
% '_',num2str(concatIndex),'.mat'));

%personal load path
load(strcat('../waveforms/UnitariesWaveforms/2013/concat_unibench_20130715_',...
    robString,'_d16/concat_unibench_20130715_',robString,'_',num2str(sequence),...
    '_',num2str(concatIndex),'.mat'));

%load('./concat_unibench_20160720_d16_qkt_fidMap_4_3.mat');

%unitary = si_calc_filt_opt_integrator_proper(opt_params,timeDelayX,timeDelayY,1);

% opt_params.rf_amp_x = 1*(1.570796326794896e+05);
% opt_params.rf_amp_y = 1*(1.570796326794896e+05);

cutoffFreqs = [0.8:0.1:5];

%filtered
[unitary,uni_hist,hist_time] = si_calc_filt_opt_integrator_proper(opt_params,cutoffFreqs(cutoff));

%unfiltered
[idealunitary,idealuni_hist,idealhist_time] = si_calc_filt_opt_integrator_proper(opt_params,1e20);

[uni_final,uniF_j,uniB_j,eigenvec_hammy_j,eigval_hammy_j]=bgrape_calc_uni_diagnostic(opt_params);
    
[outFid,uniEigFid,uniFid] = checkUniFid_fn(unitary,uni_final);

%filtered
si_unitaries = zeros(16,16,opt_params.timesteps);
transitionTimes = [1:opt_params.timesteps]*opt_params.samp_time;
for a = 1:1:opt_params.timesteps
    [compMin,compIndex] = min(abs(hist_time-transitionTimes(a)));
    si_unitaries(:,:,a) = uni_hist(:,:,compIndex);
end

%unfiltered
idealsi_unitaries = zeros(16,16,opt_params.timesteps);
idealtransitionTimes = [1:opt_params.timesteps]*opt_params.samp_time;
for a = 1:1:opt_params.timesteps
    [idealcompMin,idealcompIndex] = min(abs(hist_time-transitionTimes(a)));
    idealsi_unitaries(:,:,a) = uni_hist(:,:,idealcompIndex);
end

uniFid = zeros(opt_params.timesteps,1);
for a = 1:1:opt_params.timesteps
%     [outFid,uniEigFid,uniFid(a)] = checkUniFid_fn(si_unitaries(:,:,a),uniF_j(:,:,a));
    [outFid,uniEigFid,uniFid(a)] = checkUniFid_fn(si_unitaries(:,:,a),idealsi_unitaries(:,:,a));
end
 
si_results.int_uni = unitary;
si_results.hist_uniFid = uniFid;
si_results.hist_time = transitionTimes;
si_results.uniEigFid = uniEigFid;
% si_results.hist_uni = si_unitaries;


save(strcat('./SI_TransFuncRuns_local/si_results_TopHat_concat_unibench_20130715_',...
    robString,'_',num2str(sequence),'_',num2str(concatIndex),'_cutoffIndex_',...
    num2str(cutoff),'.mat'),'si_results');


%save('SI_TransFuncRuns/varyingKappa2kap_4step_40tot_n1_s79_0TF_si_results.mat','si_results');
%save('SI_TransFuncRuns/concat_unibench_20160720_d16_qkt_fidMap_4_3_TF_LR_RT_0Bounds_si_results.mat','si_results');
%save(strcat('./SI_TransFuncRuns/uniHist_PT_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'uni_hist');
%save(strcat('./SI_TransFuncRuns/histTime_PT_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'hist_time');
%save(strcat('./SI_TransFuncRuns/uniFid_PT_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'uniFid');
%save(strcat('./SI_TransFuncRuns/unitary_PT_',robString,'_',num2str(sequence),'_',num2str(concatIndex),'.mat'),'unitary');
%save('./ToleranceTest/uniFid.mat','uniFid');
%save('./ToleranceTest/unitary.mat','unitary');

end
