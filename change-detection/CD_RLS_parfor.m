% Change detecpplication - SOPOT and fixed-point
%
%   - An oversampled discrete white gaussian noise signal is low-passed
%   filtered with a truncated sinc function in time domain. This models a
%   generic baseband signal, such as OFDM.
%   - The considered bandwidth is 142.9 kHz.
%   - The signal is transmitted over a Rayleigh fading channel following
%   the Jakes spectrum with Doppler frequency of 104.2 Hz.
%   - A total time of 13.5 ms is considered, devided in 9 intervals with
%   the signal switched on or off at each interval.
%   At the receiver, noise is added to the faded signal.
%
%   The change detection is applied through a mean filter of the squared
%   prediction error in the adaptive filtering.
%
%   Two different scenarios are compared: (i) the SOPOT application with
%   different number of SPT terms and (ii) the fixed-point application,
%   with different word lengths.
%
%
%   Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%


clc
clear
close all

% signal definitions:
rng(5555)
ensemble = 50;  % numebr of repetitions
T = 13.5e-3;  % total time [sec]
Fs = 1e6;  % samplig rate [samp/sec]
sig_len = T*Fs;  % signal length [samp]
t = 0:(1/Fs):(T - 1/Fs);  % discrete time axis [sec]
SNR = 6.990;  % signal to noise ratio [dB]
f_d = 104.2;  % maximum Doppler shift [Hz]
filter_bw = 142.9e3;  % signal bandwidth [Hz]

% adaptive filter definitions:
N = 4;  % number of coefficients
delta = .07;  % initial value for the deterministic autocorrelation matrix
lambda = .981;  % forgetting factor
w_i = zeros(N, 1);  % initial coefficients
sptTerms = [3 4 5];  % maximum number of non-zero terms for the SOPOT
maxSums = length(sptTerms);
max_pot = 3;  % the greater power for the SOPOT
bpDepth = 16;  % maximum bit-plane depth
word_len = [7 9 11];  % word length for the fixed-point [bits]
frac_len = [2 4 6];  % depth of the representation

% generating the truncated sinc filter
t_sinc = (-7/(2*filter_bw)):1/Fs:(7/(2*filter_bw));  % sinc sample axis
trunc_sinc = (2*filter_bw)*sinc(2*filter_bw*t_sinc);  % sinc samples

% generating mask for switching on and off in time
sig_win = floor(sig_len/9);
occupied_slots_mask = zeros(1, sig_len);
occupied_slots_mask(sig_win:2*sig_win-1) = ones(1, sig_win);
occupied_slots_mask(3*sig_win:4*sig_win-1) = ones(1, sig_win);
occupied_slots_mask(5*sig_win:6*sig_win-1) = ones(1, sig_win);
occupied_slots_mask(7*sig_win:8*sig_win-1) = ones(1, sig_win);
empty_slots_mask_bool = occupied_slots_mask == 0;  % true when empty
occupied_slots_mask_bool = ~empty_slots_mask_bool;
total_empty_time = 5*sig_win;
total_occupied_time = 4*sig_win;

% for the ROC curve
thres = 0:.00005:.2;  % threshold axis
H = length(thres);  % for threshold index

% setup for fixed-point arithmetic
G = globalfimath;
G.ProductMode = 'SpecifyPrecision';
G.SumMode = 'SpecifyPrecision';

% setup structs for the operations
S_sopot = struct('filterOrderNo', N-1, 'initialCoefficients', w_i,...
    'lambda', lambda, 'delta', delta, 'windowSize', 30, 'bitplaneDepth',...
    bpDepth);
S_fxdpt = struct('filterOrderNo', N-1, 'initialCoefficients', w_i,...
                 'lambda', lambda, 'delta', delta, 'windowSize', 30);

% memory allocation
sopot_alarm_out = zeros(H, sig_len, maxSums, ensemble);
fxdpt_alarm_out = zeros(H, sig_len, maxSums, ensemble);
CD_time_sopot = zeros(H, 8, maxSums, ensemble);
CD_time_fxdpt = zeros(H, 8, maxSums, ensemble);
fp_ratio_sopot = zeros(H, maxSums, ensemble);
fp_ratio_fxdpt = zeros(H, maxSums, ensemble);
fn_ratio_sopot = zeros(H, maxSums, ensemble);
fn_ratio_fxdpt = zeros(H, maxSums, ensemble);
sopot_rate_fp = zeros(H, maxSums, ensemble);
sopot_rate_fn = zeros(H, maxSums, ensemble);
fxdpt_rate_fp = zeros(H, maxSums, ensemble);
fxdpt_rate_fn = zeros(H, maxSums, ensemble);

temp_count = zeros(maxSums, 1);

% computing
for m = 1:ensemble  % monte carlo loop
%     tic;
    % generating the limited-band signal
    rand_sig = randn(1, sig_len);
    filt_sig = filter(trunc_sinc, 1, rand_sig);
    filt_sig = filt_sig/norm(filt_sig, inf);
    % switch on and off in time
    trans_sig = filt_sig.*occupied_slots_mask;  % signal to be transmitted
    
    % generating the Rayleigh fading channel
    raychan = comm.RayleighChannel('SampleRate', Fs,...
        'MaximumDopplerShift', f_d);
    
    % transmitting through the Rayleigh channel
    fin_sig = real(raychan(trans_sig').');
    
    % adding noise
    P_sig = sum(abs(fin_sig).^2, 2)./sig_len;  % signal power
    P_noise = P_sig.*(10.^(-.1*SNR));  % new noise power
    noise = randn(1, sig_len);  % noise generation
    dummy = sum(noise.^2, 2)./sig_len;  %  old noise power
    noise = noise.*sqrt(P_noise/dummy);  % power adjustment
    reciv_sig = fin_sig + noise;  % signal at the reciever
    
    for spt_iter = 1:maxSums  % the number of active bits loop
        % adjusting the global fixed-point arithmetic
        G = globalfimath('ProductWordLength', word_len(spt_iter),...
            'ProductFractionLength', frac_len(spt_iter),...
            'SumWordLength', word_len(spt_iter),...
            'SumFractionLength', frac_len(spt_iter));
        temp = sptTerms(spt_iter); 
        wordLength = word_len(spt_iter);
        fracLength = frac_len(spt_iter);
        parfor h = 1:H  % the different thresholds loop
            tic
            % broadcasted variables:
            wordLength_pf = wordLength;
            fracLength_pf = fracLength;
            sigIn = reciv_sig';
            test(h) = thres(h);
            % adaptive filtering
            % SOPOT
            sopot_alarm(h, :) = cd_rls_sopot(sigIn, temp, max_pot,...
                                       temp, thres(h), S_sopot);
            test2(h) = thres(h);
            % fixed-point
            fxdpt_alarm(h, :) = cd_rls_fxdpt(sigIn, wordLength_pf,...
                                        fracLength_pf, thres(h),...
                                        S_fxdpt);
            toc
        end
        
%         sopot_alarm_out(:, :, spt_iter, m) = sopot_alarm;
%         fxdpt_alarm_out(:, :, spt_iter, m) = fxdpt_alarm;
        
        for it = 1:H
            [CD_time_sopot(it, :, spt_iter, m),...
                fp_ratio_sopot(it, spt_iter, m),...
                fn_ratio_sopot(it, spt_iter, m)] = alarm_proc(...
                sopot_alarm(it, :), empty_slots_mask_bool, 30);
            [CD_time_fxdpt(it, :, spt_iter, m),...
                fp_ratio_fxdpt(it, spt_iter, m),...
                fn_ratio_fxdpt(it, spt_iter, m)] = alarm_proc(...
                fxdpt_alarm(it, :), empty_slots_mask_bool, 30);
        end
        sav_temp_count(spt_iter) = temp;
        % No need for this now.
        sopot_rate_fp_aux = zeros(H, 1);
        sopot_rate_fn_aux = zeros(H, 1);
        fxdpt_rate_fp_aux = zeros(H, 1);
        fxdpt_rate_fn_aux = zeros(H, 1);
        
        for it = 1:H
        % classify data
        % positive: empty channel
        % negative: occuped channel
            
        sopot_dif = xor(sopot_alarm(it, :), empty_slots_mask_bool);
        fxdpt_dif = xor(fxdpt_alarm(it, :), empty_slots_mask_bool);
        
        % FALSE POSITIVE
        sopot_fp = sopot_dif(empty_slots_mask_bool);
        sopot_rate_fp_aux(it) = sum(double(sopot_fp), 2)./...
            total_empty_time;
        fxdpt_fp = fxdpt_dif(empty_slots_mask_bool);
        fxdpt_rate_fp_aux(it) = sum(double(fxdpt_fp), 2)./...
            total_empty_time;
        
        % FALSE NEGATIVE
        sopot_fn = sopot_dif(occupied_slots_mask_bool);
        sopot_rate_fn_aux(it) = sum(double(sopot_fn), 2)./...
            total_occupied_time;
        fxdpt_fn = fxdpt_dif(occupied_slots_mask_bool);
        fxdpt_rate_fn_aux(it) = sum(double(fxdpt_fn), 2)./...
            total_occupied_time;
        end
        
        sopot_rate_fp(:, spt_iter, m) = sopot_rate_fp_aux;
        sopot_rate_fn(:, spt_iter, m) = sopot_rate_fn_aux;
        fxdpt_rate_fp(:, spt_iter, m) = fxdpt_rate_fp_aux;
        fxdpt_rate_fn(:, spt_iter, m) = fxdpt_rate_fn_aux;
%         
    end
%     if m == 1
%         time = toc;
%         sec = floor((time./60 - floor(time./60)).*60);
%         min_aux = floor(time./60);
%         min = (min_aux./60 - floor(min_aux./60)).*60;
%         hour_aux = floor(min_aux./60);
%         hour = (hour_aux./24 - floor(hour_aux./24)).*24;
%         day = floor(hour_aux./24);
%         time_est = toc.*ensemble;
%         sec_est = floor((time_est./60 - floor(time_est./60)).*60);
%         min_aux = floor(time_est./60);
%         min_est = (min_aux./60 - floor(min_aux./60)).*60;
%         hour_aux = floor(min_aux./60);
%         hour_est = (hour_aux./24 - floor(hour_aux./24)).*24;
%         day_est = floor(hour_aux./24);
%     else
%         time = time + toc;
%         sec = floor((time./60 - floor(time./60)).*60);
%         min_aux = floor(time./60);
%         min = (min_aux./60 - floor(min_aux./60)).*60;
%         hour_aux = floor(min_aux./60);
%         hour = (hour_aux./24 - floor(hour_aux./24)).*24;
%         day = floor(hour_aux./24);
%     end
%     fprintf('O tempo total decorrido eh: %u dias, %u horas, %u minutos e %u segundos.\n',...
%         day, hour, min, sec)
%     fprintf('O tempo total estimado eh: %u dias, %u horas, %u minutos e %u segundos.\n',...
%         day_est, hour_est, min_est, sec_est);
end

sav_CD_sopot = CD_time_sopot;
sav_CD_fxdpt = CD_time_fxdpt;
sav_fp_ratio_sopot_mean = sum(fp_ratio_sopot, 3)./ensemble;
sav_fp_ratio_fxdpt_mean = sum(fp_ratio_fxdpt, 3)./ensemble;
sav_fn_ratio_sopot_mean = sum(fn_ratio_sopot, 3)./ensemble;
sav_fn_ratio_fxdpt_mean = sum(fn_ratio_fxdpt, 3)./ensemble;


% % the means
rate_sopot_fp_mean = sum(sopot_rate_fp, 3)./ensemble;
rate_sopot_fn_mean = sum(sopot_rate_fn, 3)./ensemble;
rate_fxdpt_fp_mean = sum(fxdpt_rate_fp, 3)./ensemble;
rate_fxdpt_fn_mean = sum(fxdpt_rate_fn, 3)./ensemble;
% 
% % saving variables
save('results_CD_RLS.mat', 'rate*');

% save('CD_RLS_alarms.mat', 'sopot_alarm_out', 'fxdpt_alarm_out')
save('results_test.mat', 'sav_*')

% EoF


