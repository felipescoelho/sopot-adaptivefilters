% SOPOT project, offline application for LMS, NLMS and RLS
% SYSTEM IDENTIFICATION
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%

clc
clear
close all

% Adjustables:
Adders_max = 32;  % Max number of bits in the approximation
N = 4;  % Max number of coefficients in the filters
H_dummy = randn(N, 1);
H = H_dummy./norm(H_dummy);  % Unknown system
lambda = .97;  % Forgeting factor
delta = .42;  % small const. to initialize the Hermitian
mu_LMS = .1;  % step size for the LMS filter
mu_NLMS = .098;  % step size for the NLMS filter
gamma = 1e-12;  % Small constant to avoid singularity

K = 1000;  % No. iterations
ensemble = 100;  % Monte Carlo ensemble

save('meta.mat', 'lambda', 'delta', 'mu_*', 'gamma', 'K', 'Adders_max')

% Memory allocation:


% Computing:
P = floor(sqrt(N));  % P parameter for the MPGBP algorithm
parfor iter = 1:ensemble
    % Creating signals
    x = randn(K, 1);  % Input signal
    sigma_x2 = var(x);  % Signal power
    x = x./sqrt(sigma_x2);  % Signal power normalization
    
    d_aux = filter(H, 1, x);  % Power adjusted interference signal
    Pd_aux = (norm(d_aux, 2).^2)./K;  % d_aux power
    n = randn(K, 1);  % corruptive white noise
    sigma_n2 = var(n);
    
    SNR = [-3 3 10 15];  % SNR between n and x [dB]
    snr_len = length(SNR);  % length of the SNR axis
    
    % Memory allocation:
    km_RLS = zeros(Adders_max, snr_len);
    km_LMS = zeros(Adders_max, snr_len);
    km_NLMS = zeros(Adders_max, snr_len);
    w_RLS = zeros(N, K+1, snr_len);
    w_RLS_SPT = zeros(N, K+1, Adders_max, snr_len);
    w_LMS = zeros(N, K+1, snr_len);
    w_LMS_SPT = zeros(N, K+1, Adders_max, snr_len);
    w_NLMS = zeros(N, K+1, snr_len);
    w_NLMS_SPT = zeros(N, K+1, Adders_max, snr_len);
    norm_w_RLS = zeros(K+1, snr_len);
    norm_w_RLS_SPT = zeros(K+1, Adders_max, snr_len);
    norm_w_LMS = zeros(K+1, snr_len);
    norm_w_LMS_SPT = zeros(K+1, Adders_max, snr_len);
    norm_w_NLMS = zeros(K+1, snr_len);
    norm_w_NLMS_SPT = zeros(K+1, Adders_max, snr_len);
    MSD_RLS = zeros(K+1, snr_len);
    MSD_RLS_SPT = zeros(K+1, Adders_max, snr_len);
    MSE_RLS = zeros(K, snr_len);
    MSE_RLS_SPT = zeros(K, Adders_max, snr_len);
    MSD_LMS = zeros(K+1, snr_len);
    MSD_LMS_SPT = zeros(K+1, Adders_max, snr_len);
    MSE_LMS = zeros(K, snr_len);
    MSE_LMS_SPT = zeros(K, Adders_max, snr_len);
    MSD_NLMS = zeros(K+1, snr_len);
    MSD_NLMS_SPT = zeros(K+1, Adders_max, snr_len);
    MSE_NLMS = zeros(K, snr_len);
    MSE_NLMS_SPT = zeros(K, Adders_max, snr_len);
    
    % Computing
    for snr = 1:snr_len
        W = zeros(N, 1);  % Intial coefficients vectors
        Pn = Pd_aux.*(10.^(-.1.*SNR(snr)));  % Noise power
        n = sqrt(Pn./sigma_n2).*n;  % Noise power adjustment
        d = d_aux + n;  % Desired signal
        
        % Setting up structs:
        S_RLS = struct('filterOrderNo', (N-1), 'delta', delta, 'lambda', lambda);
        S_LMS = struct('filterOrderNo', (N-1), 'step', mu_LMS, 'initialCoefficients', W);
        S_NLMS = struct('step', mu_NLMS, 'filterOrderNo', (N-1), 'initialCoefficients', W, 'gamma', gamma);
        
        % Adaptive filtering process:
        [~, e_RLS, w_RLS(:, :, snr), ~, ~] = RLS(d, transpose(x), S_RLS);
        [~, e_LMS, w_LMS(:, :, snr)] = LMS(d, transpose(x), S_LMS);
        [~, e_NLMS, w_NLMS(:, :, snr)] = NLMS(d, transpose(x), S_NLMS);
        
        for nAdders = 1:Adders_max
            
            % Setting up structs, for SOPOT:
            S_RLS_SPT = struct('filterOrderNo', (N-1), 'delta', delta,...
                'lambda', lambda, 'adders', nAdders);
            S_LMS_SPT = struct('step', mu_LMS, 'filterOrderNo', (N-1),...
                'initialCoefficients', W, 'adders', nAdders);
            S_NLMS_SPT = struct('step', mu_NLMS, 'filterOrderNo', (N-1), 'initialCoefficients', W, 'gamma', gamma, 'adders', nAdders);
            
            % Adaptive filtering process, for SOPOT:
            [~, e_RLS_SPT, w_RLS_SPT(:, :, nAdders, snr), ~, ~, km_RLS(nAdders, snr)] = RLS_SOPOT(d, transpose(x), S_RLS_SPT);
            [~, e_LMS_SPT, w_LMS_SPT(:, :, nAdders, snr), km_LMS(nAdders, snr)] = LMS_SOPOT(d, transpose(x), S_LMS_SPT);
            [~, e_NLMS_SPT, w_NLMS_SPT(:, :, nAdders, snr), km_NLMS(nAdders, snr)] = NLMS_SOPOT(d, transpose(x), S_NLMS_SPT);
            
            % MSE and MSD
            for k = 1:K+1
                norm_w_RLS_SPT(k, nAdders, snr) = norm(w_RLS_SPT(:, k, nAdders, snr));
                norm_w_LMS_SPT(k, nAdders, snr) = norm(w_LMS_SPT(:, k, nAdders, snr));
                norm_w_NLMS_SPT(k, nAdders, snr) = norm(w_NLMS_SPT(:, k, nAdders, snr));
            end
            deviation_RLS_SPT = 1 - norm_w_RLS_SPT(:, nAdders, snr)
            MSD_RLS_SPT(:, nAdders, snr) = MSD_RLS_SPT(:, nAdders, snr) + abs(deviation_RLS_SPT).^2;
            MSE_RLS_SPT(:, nAdders, snr) = MSE_RLS_SPT(:, nAdders, snr) + abs(e_RLS_SPT).^2;
            deviation_LMS_SPT = 1 - norm_w_LMS_SPT(:, nAdders, snr)
            MSD_LMS_SPT(:, nAdders, snr) = MSD_LMS_SPT(:, nAdders, snr) + abs(deviation_LMS_SPT).^2;
            MSE_LMS_SPT(:, nAdders, snr) = MSE_LMS_SPT(:, nAdders, snr) + abs(e_LMS_SPT).^2;
            deviation_NLMS_SPT = 1 - norm_w_NLMS_SPT(:, nAdders, snr)
            MSD_NLMS_SPT(:, nAdders, snr) = MSD_NLMS_SPT(:, nAdders, snr) + abs(deviation_NLMS_SPT).^2;
            MSE_NLMS_SPT(:, nAdders, snr) = MSE_NLMS_SPT(:, nAdders, snr) + abs(e_NLMS_SPT).^2;
        end
        for k = 1:K+1
            norm_w_RLS(k, snr) = norm(w_RLS(:, k, snr));
            norm_w_LMS(k, snr) = norm(w_LMS(:, k, snr));
            norm_w_NLMS(k, snr) = norm(w_NLMS(:, k, snr));
        end
        deviation_RLS = 1 - norm_w_RLS(:, snr)
        MSD_RLS(:, snr) = MSD_RLS(:, snr) + abs(deviation_RLS).^2;
        MSE_RLS(:, snr) = MSE_RLS(:, snr) + abs(e_RLS).^2;
        deviation_LMS = 1 - norm_w_LMS(:, snr)
        MSD_LMS(:, snr) = MSD_LMS(:, snr) + abs(deviation_LMS).^2;
        MSE_LMS(:, snr) = MSE_LMS(:, snr) + abs(e_LMS).^2;
        deviation_NLMS = 1 - norm_w_NLMS(:, snr)
        MSD_NLMS(:, snr) = MSD_NLMS(:, snr) + abs(deviation_NLMS).^2;
        MSE_NLMS(:, snr) = MSE_NLMS(:, snr) + abs(e_NLMS).^2;
    end
    km_RLS_out(:, :, iter) = km_RLS;
    km_LMS_out(:, :, iter) = km_LMS;
    km_NLMS_out(:, :, iter) = km_NLMS;
    w_RLS_out(:, :, :, iter) = w_RLS;
    w_RLS_SPT_out(:, :, :, :, iter) = w_RLS_SPT;
    w_LMS_out(:, :, :, iter) = w_LMS;
    w_LMS_SPT_out(:, :, :, :, iter) = w_LMS_SPT;
    w_NLMS_out(:, :, :, iter) = w_NLMS;
    w_NLMS_SPT_out(:, :, :, :, iter) = w_NLMS_SPT;
    norm_w_RLS_out(:, :, iter) = norm_w_RLS;
    norm_w_RLS_SPT_out(:, :, :, iter) = norm_w_RLS_SPT;
    norm_w_LMS_out(:, :, iter) = norm_w_LMS;
    norm_w_LMS_SPT_out(:, :, :, iter) = norm_w_LMS_SPT;
    norm_w_NLMS_out(:, :, iter) = norm_w_NLMS;
    norm_w_NLMS_SPT_out(:, :, :, iter) = norm_w_NLMS_SPT;
    MSD_RLS_out(:, :, iter) = MSD_RLS;
    MSD_RLS_SPT_out(:, :, :, iter) = MSD_RLS_SPT;
    MSE_RLS_out(:, :, iter) = MSE_RLS;
    MSE_RLS_SPT_out(:, :, :, iter) = MSE_RLS_SPT;
    MSD_LMS_out(:, :, iter) = MSD_LMS;
    MSD_LMS_SPT_out(:, :, :, iter) = MSD_LMS_SPT;
    MSE_LMS_out(:, :, iter) = MSE_LMS;
    MSE_LMS_SPT_out(:, :, :, iter) = MSE_LMS_SPT;
    MSD_NLMS_out(:, :, iter) = MSD_NLMS;
    MSD_NLMS_SPT_out(:, :, :, iter) = MSD_NLMS_SPT;
    MSE_NLMS_out(:, :, iter) = MSE_NLMS;
    MSE_NLMS_SPT_out(:, :, :, iter) = MSE_NLMS_SPT;
    
end

% Average calc.:
km_RLS_avg = sum(km_RLS_out, 3)./ensemble;
km_LMS_avg = sum(km_LMS_out, 3)./ensemble;
km_NLMS_avg = sum(km_NLMS_out, 3)./ensemble;
w_RLS_avg = sum(w_RLS_out, 4)./ensemble;
w_RLS_SPT_avg = sum(w_RLS_SPT_out, 5)./ensemble;
w_LMS_avg = sum(w_LMS_out, 4)./ensemble;
w_LMS_SPT_avg = sum(w_LMS_SPT_out, 5)./ensemble;
w_NLMS_avg = sum(w_NLMS_out, 4)./ensemble;
w_NLMS_SPT_avg = sum(w_NLMS_SPT_out, 5)./ensemble;
norm_w_RLS_avg = sum(norm_w_RLS_out, 3)./ensemble;
norm_w_RLS_SPT_avg = sum(norm_w_RLS_SPT_out, 4)./ensemble;
norm_w_LMS_avg = sum(norm_w_LMS_out, 3)./ensemble;
norm_w_LMS_SPT_avg = sum(norm_w_LMS_SPT_out, 4)./ensemble;
norm_w_NLMS_avg = sum(norm_w_NLMS_out, 3)./ensemble;
norm_w_NLMS_SPT_avg = sum(norm_w_NLMS_SPT_out, 4)./ensemble;
MSD_RLS_avg = sum(MSD_RLS_out, 3)./ensemble;
MSD_RLS_SPT_avg = sum(MSD_RLS_SPT_out, 4)./ensemble;
MSE_RLS_avg = sum(MSE_RLS_out, 3)./ensemble;
MSE_RLS_SPT_avg = sum(MSE_RLS_SPT_out, 4)./ensemble;
MSD_LMS_avg = sum(MSD_LMS_out, 3)./ensemble;
MSD_LMS_SPT_avg = sum(MSD_LMS_SPT_out, 4)./ensemble;
MSE_LMS_avg = sum(MSE_LMS_out, 3)./ensemble;
MSE_LMS_SPT_avg = sum(MSE_LMS_SPT_out, 4)./ensemble;
MSD_NLMS_avg = sum(MSD_NLMS_out, 3)./ensemble;
MSD_NLMS_SPT_avg = sum(MSD_NLMS_SPT_out, 4)./ensemble;
MSE_NLMS_avg = sum(MSE_NLMS_out, 3)./ensemble;
MSE_NLMS_SPT_avg = sum(MSE_NLMS_SPT_out, 4)./ensemble;


save('systemID_off.mat', '*_avg')


% EoF

