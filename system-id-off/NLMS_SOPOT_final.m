    %                NLMS/SOPOT - Off-line Filter Approximation and In-loop


clc
close all
clear all                        %#ok<CLALL>


%    File Names:
original          = '../Originais/handel.wav';
ruido             = '../Originais/n1.wav';
corrompido        = '../Processados/Corrompido.wav';
filtrado_NLMS     = '../Processados/limpo_NLMS.wav';
filtrado_SOPOT    = '../Processados/limpo_SOPOT.wav';


%    Adjustables:
Adders_max          = 22;           % Max number of Bits in the Approximation
N_max               = 5;            % Max number of coefficients ion the Filters
mu                  = 0.0019;       % NLMS convergence step size
gamma               = 1e-12;        % small positive interger to avoid singularity
SNR_n1              = 10;           % SNR between n1 and x
SNR_n2              = 10;           % SNR between n2 and x


%    Reading the orignal audio file:
[x, Fs]  = audioread(original);             % x is the reference signal and Fs is its sample rate

[~, ncolumn] = size(x);
if ncolumn == 2
    x    = x(:, 1);
end

L        = length(x);                       % reference signal distribution length




%%%%%%%%%%%%%%%%%    Reducao do numero de iteracoes para ir mais rapido.
L = floor((3*L)/7);
x = x(1:L);





%    Generating the interference noise n1:
rng('default');
n1 = randn(L, 1);     % interference signal


%    Adding the interference noise to the reference signal:
Px         = (norm(x, 2) ^ 2) / L;
Pn1        = Px * (10 ^ (- 0.1 * SNR_n1));
Pn1_aux    = (norm(n1, 2) ^ 2) / L;
n1         = sqrt(Pn1 / Pn1_aux) * n1;
corromp    = n1 + x;
aux        = norm(corromp, inf);

if aux > 1
    corromp    = corromp / aux;
    x2         = x / aux;
    n1         = n1 / aux;
else
    x2 = x;
end

Px2       = (norm(x2, 2)^2)/L;
audiowrite(corrompido, corromp, Fs);


%    Adjusting the input noise n2:
a          = [1 -0.95];                             % coeficientes do denominador da funcao 'filter'
b          = 1;                                     % coeficientes do numerador da funcao 'filter'
n2         = filter(b, a, n1);                          
Pn2        = Px * (10 ^ (- 0.1 * SNR_n2));
Pn2_aux    = (norm(n2, 2) ^ 2) / L;
n2         = sqrt(Pn2 / Pn2_aux) * n2;


%    Initializing & Allocating memory:
wn                  = zeros(N_max-1, N_max);                     % coefficient vector for each iteration and realization
d                   = corromp;                                   % desired signal
ho                  = zeros(Adders_max, N_max);
hn                  = zeros(Adders_max, N_max, N_max-1);
w_spt_addn          = zeros(Adders_max, N_max, N_max-1);
w_spt_add           = zeros(Adders_max, N_max);
km2                 = zeros(Adders_max, 1);
km_inloop           = zeros(Adders_max, 1);
km_inloopn          = zeros(N_max-1, Adders_max);
km_n                = zeros(N_max-1, Adders_max);
ruid_sopot          = zeros(L, Adders_max, N_max-1);
e_spt_aux           = zeros(L, Adders_max, N_max-1);
ruid_spt_inloop     = zeros(L, Adders_max, N_max-1);
e_NLMS_aux          = zeros(L, Adders_max, N_max-1);
ruid_NLMS_offline   = zeros(L, Adders_max, N_max-1);
w_coeff_add_N       = zeros(N_max, L+1, Adders_max, N_max-1);
w_spt_coeff_add_N   = zeros(N_max, L+1, Adders_max, N_max-1);


%    Computing:
for N = 2:N_max
      
    P        = floor(sqrt(N));
    W        = zeros(N, (L+1));
       
    
    for nAdders = 1:Adders_max
    
        S1 = struct('step', mu, 'filterOrderNo', (N-1), 'initialCoefficients', W(:, 1), 'gamma', gamma, 'adders', nAdders);
        
        S2 = struct('step', mu, 'filterOrderNo', (N-1), 'initialCoefficients', W(:, 1), 'gamma', gamma);
        
        [y, e, w]                               = NLMS(d, transpose(n2), S2);
        
        [y_spt, e_spt, w_spt, km]               = NLMS_SOPOTinloop(d, transpose(n2), S1);
        
        km_inloop(nAdders, 1)                   = km(end);
        
        w_sptAUX                                = w_spt(:, end);
        
        wo                                      = w(:, end);
        
        [h_aux, km2(nAdders, 1)]                = mpgbp_filter_approxP(wo, nAdders, P);
        
        ruid_sopot(:, nAdders, N-1)             = n1 - filter(h_aux, 1, n2);
        
        e_spt_aux(:, nAdders, N-1)              = abs(e_spt.^2);
        
        ruid_spt_inloop(:, nAdders, N-1)        = n1 - filter(w_sptAUX, 1, n2);
        
        e_NLMS_aux(:, nAdders, N-1)             = abs(e.^2);
        
        ruid_NLMS_offline(:, nAdders, N-1)      = n1 - filter(wo, 1, n2);
        
        
        
        for it = 1:N
            
            w_coeff_add_N(it, :, nAdders, N-1)      = w(it, :);
            
            w_spt_add(nAdders, it)                  = w_sptAUX(it);
            
            ho(nAdders, it)                         = h_aux(it);
            
            w_spt_coeff_add_N(it, :, nAdders, N-1)  = w_spt(it, :);
            
            wn(N-1, it)                             = wo(it);
            
        end
    end
    
    km_n (N-1, :)         = km2;
    
    km_inloopn(N-1, :)    = km_inloop;
    
    w_spt_addn(:, :, N-1) = w_spt_add;
    
    hn(:, :, N-1)         = ho;
    
end



%    Avaliacao de Desempenho da aproximacao por SOPOT MSE e Analise da anulacao do efeito da AR de primeira Ordem:
Wn                                  = zeros(N_max-1, 512);
Hn                                  = zeros(N_max-1, 512, Adders_max);
Wspt                                = zeros(N_max-1, 512, Adders_max);
Resposta_cascata_NLMS               = zeros(N_max-1, 512);
Resposta_cascata_SPT                = zeros(Adders_max, 512, N_max-1);
Resposta_cascata_SPTinloop          = zeros(Adders_max, 512, N_max-1);
MSE_aux_NLMSandApprox               = zeros(N_max-1, Adders_max);
MSE_aux_cascade_NLMSandApprox       = zeros(N_max-1, Adders_max);
MSE_aux_cascade_NLMSandInloop       = zeros(N_max-1, Adders_max);
MSE_aux_cascade_ApproxandInloop     = zeros(N_max-1, Adders_max);
MSE_aux_NLMSandInloop               = zeros(N_max-1, Adders_max);
MSE_aux_ApproxandInloop             = zeros(N_max-1, Adders_max);
Adders_perCoeff                     = zeros(N_max-1, Adders_max);
tabela_aux_adders                   = zeros(N_max-1, Adders_max);
tabela_aux_coeff                    = zeros(N_max-1, Adders_max);
SNR_SOPOT_aux                       = zeros(N_max-1, Adders_max);
SNR_NLMS_offline_aux                = zeros(N_max-1, Adders_max);
SNR_NLMS_loop_aux                   = zeros(N_max-1, Adders_max);
SNR_spt_Inloop_aux                  = zeros(N_max-1, Adders_max);


for N = 1:N_max-1
        
    Wn(N, :)                    = freqz(wn(N, :), 1);
    
    Resposta_cascata_NLMS(N, :) = freqz(wn(N, :), a);
        
    for nAdders = 1:Adders_max
            
        Hn(N, :, nAdders)                           = freqz(hn(nAdders, :, N), 1);
        
        Wspt(N, :, nAdders)                         = freqz(w_spt_addn(nAdders, :, N), 1);
            
        Resposta_cascata_SPT(N, :, nAdders)         = freqz(hn(nAdders, :, N), a);
        
        Resposta_cascata_SPTinloop(N, :, nAdders)   = freqz(w_spt_addn(nAdders, :, N), a);
            
        MSE_aux_NLMSandApprox(N, nAdders)           = (sum(((Wn(N, :) - Hn(N, :, nAdders)).^2), 2))/512;
        
        MSE_aux_NLMSandInloop(N, nAdders)           = (sum(((Wn(N, :) - Wspt(N, :, nAdders)).^2), 2))/512;
        
        MSE_aux_ApproxandInloop(N, nAdders)         = (sum(((Hn(N, :, nAdders) - Wspt(N, :, nAdders)).^2), 2))/512;
        
        MSE_aux_cascade_NLMSandApprox(N, nAdders)   = (sum(((Resposta_cascata_NLMS(N, :) - Resposta_cascata_SPT(N, :, nAdders)).^2), 2))/512;
        
        MSE_aux_cascade_NLMSandInloop(N, nAdders)   = (sum(((Resposta_cascata_NLMS(N, :) - Resposta_cascata_SPTinloop(N, :, nAdders)).^2), 2))/512;
        
        MSE_aux_cascade_ApproxandInloop(N, nAdders) = (sum(((Resposta_cascata_SPT(N, :, nAdders) - Resposta_cascata_SPTinloop(N, :, nAdders)).^2), 2))/512;
            
        Adders_perCoeff(N, nAdders)                 = nAdders/(N+1);
            
        tabela_aux_adders(N, nAdders)               = nAdders;
            
        tabela_aux_coeff(N, nAdders)                = N+1;
        
        SNR_SOPOT_aux(N, nAdders)                   = Px2/((norm(ruid_sopot(:, nAdders, N), 2)^2)/L);
        
        SNR_NLMS_offline_aux(N, nAdders)            = Px2/((norm(ruid_NLMS_offline(:, nAdders, N), 2)^2)/L);
        
        SNR_spt_Inloop_aux(N, nAdders)              = Px2/((norm(ruid_spt_inloop(:, nAdders, N), 2)^2)/L);
            
    end
end


g                           = (N_max-1) * Adders_max;
SNR_NLMS_offline            = 10*log10(SNR_NLMS_offline_aux);
SNR_SOPOT                   = 10*log10(SNR_SOPOT_aux);
SNR_spt_Inloop              = 10*log10(SNR_spt_Inloop_aux);
MSE_NLMSandApprox           = 20*log10(abs(MSE_aux_NLMSandApprox));
MSE_cascade_NLMSandApprox   = 20*log10(abs(MSE_aux_cascade_NLMSandApprox));
MSE_cascade_NLMSandInloop   = 20*log10(abs(MSE_aux_cascade_NLMSandInloop));
MSE_cascade_ApproxandInloop = 20*log10(abs(MSE_aux_cascade_ApproxandInloop));
MSE_NLMSandInloop           = 20*log10(abs(MSE_aux_NLMSandInloop));
MSE_ApproxandInloop         = 20*log10(abs(MSE_aux_ApproxandInloop));

tabela_aux3(:, 1)   = tabela_aux_adders(:);
tabela_aux3(:, 2)   = tabela_aux_coeff(:);
tabela_aux3(:, 3)   = Adders_perCoeff(:);
tabela_aux3(:, 4)   = -km_n(:);
tabela_aux3(:, 5)   = -km_inloopn(:);
tabela_aux3(:, 6)   = MSE_NLMSandApprox(:);
tabela_aux3(:, 7)   = MSE_cascade_NLMSandApprox(:);
tabela_aux3(:, 8)   = MSE_cascade_NLMSandInloop(:);
tabela_aux3(:, 9)   = MSE_cascade_ApproxandInloop(:);
tabela_aux3(:, 10)  = MSE_NLMSandInloop(:);
tabela_aux3(:, 11)  = MSE_ApproxandInloop(:);
tabela_aux3(:, 12)  = SNR_NLMS_offline(:);
tabela_aux3(:, 13)  = SNR_SOPOT(:);
tabela_aux3(:, 14)  = SNR_spt_Inloop(:);

tabela_aux4 = sortrows(tabela_aux3);

for N = 1:g
    
    if tabela_aux4(N,1) == 1
        
        break
        
    end
end

tabela = zeros(g-N, 14);

for nAdders = N:g
    
    tabela(nAdders-N+1,:) = tabela_aux4(nAdders,:);
    
end

%%%%%%%%% Grafico km, off-line/inloop:

N = [1 2 3 4];

tabela_km = sortrows(tabela_aux3, 2);
curva_km = zeros(Adders_max, 3, 4);

for n = 1:4
    
    curva_km(:, 1, n) = tabela_km(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 1);
    curva_km(:, 2, n) = tabela_km(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 4);
    curva_km(:, 3, n) = tabela_km(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 5);
    
end

figure(1),...
    p1 = plot(curva_km(:,1,1), curva_km(:,2,1)); hold on,...
    p2 = plot(curva_km(:,1,2), curva_km(:,2,2));...
    p3 = plot(curva_km(:,1,3), curva_km(:,2,3));...
    p4 = plot(curva_km(:,1,4), curva_km(:,2,4));...
    ax = gca;
    ax.ColorOrderIndex = 1;
    p5 = plot(curva_km(:,1,1), curva_km(:,3,1), '--');...
    p6 = plot(curva_km(:,1,2), curva_km(:,3,2), '--');...
    p7 = plot(curva_km(:,1,3), curva_km(:,3,3), '--');...
    p8 = plot(curva_km(:,1,4), curva_km(:,3,4), '--');...
    legend('Off-line - N = 2', 'Off-line - N = 3', 'Off-line - N = 4', 'Off-line - N = 5', 'In-loop - N = 2', 'In-loop - N = 3', 'In-loop - N = 4', 'In-loop - N = 5'),...
    xlabel('Adders'), ylabel('Bit-depth, km'),...
    grid on, grid minor, hold off;
    p1.LineWidth = 1;
    p2.LineWidth = 1;
    p3.LineWidth = 1;
    p4.LineWidth = 1;
    p5.LineWidth = 1;
    p6.LineWidth = 1;
    p7.LineWidth = 1;
    p8.LineWidth = 1;
    

%%%%%%%%%%%%%%%%%%%%  Graficos (1) MSE - Adders:

%  MSE - NLMS/SOPOT e NLMS/Inloop-SOPOT

N  = [1 2 3 4];

curva1          = zeros(Adders_max, 2, 4);
curva2          = zeros(Adders_max, 2, 4);
curva3          = zeros(Adders_max, 2, 4);
curva4          = zeros(Adders_max, 2, 4);
curva_aux1      = zeros(Adders_max, 2, 4);
curva_aux2      = zeros(Adders_max, 2, 4);
curva_aux3      = zeros(Adders_max, 2, 4);
curva_aux4      = zeros(Adders_max, 2, 4);
nnzer1          = zeros(1, 4);
nnzer2          = zeros(1, 4);
nnzer3          = zeros(1, 4);
nnzer4          = zeros(1, 4);

tabela_grafic   = sortrows(tabela, 2);

for nAdders = 1:4
    
    curva1(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 1);    
    curva1(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 6);
    
    curva2(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 1);    
    curva2(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 10);
    
    curva3(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 1);    
    curva3(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 7);
    
    curva4(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 1);    
    curva4(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 8);
    

    for k = 2:Adders_max
        if curva1(k, 2, nAdders) >= curva1(k-1, 2, nAdders)
        else
            curva_aux1(k, :, nAdders) = curva1(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva2(k, 2, nAdders) >= curva2(k-1, 2, nAdders)
        else
            curva_aux2(k, :, nAdders) = curva2(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva3(k, 2, nAdders) >= curva3(k-1, 2, nAdders)
        else
            curva_aux3(k, :, nAdders) = curva3(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva4(k, 2, nAdders) >= curva4(k-1, 2, nAdders)
        else
            curva_aux4(k, :, nAdders) = curva4(k, :, nAdders);
        end
    end   
    
    nnzer1(nAdders)              =   nnz(curva_aux1(:, 1, nAdders));
    nnzer2(nAdders)              =   nnz(curva_aux2(:, 1, nAdders));
    nnzer3(nAdders)              =   nnz(curva_aux3(:, 1, nAdders));
    nnzer4(nAdders)              =   nnz(curva_aux4(:, 1, nAdders));
    curva_aux1(:, :, nAdders)    =   sortrows(curva_aux1(:, :, nAdders), 1);
    curva_aux2(:, :, nAdders)    =   sortrows(curva_aux2(:, :, nAdders), 1);
    curva_aux3(:, :, nAdders)    =   sortrows(curva_aux3(:, :, nAdders), 1);
    curva_aux4(:, :, nAdders)    =   sortrows(curva_aux4(:, :, nAdders), 1);
    

end

curv_1   =  zeros(nnzer1(1), 2);
curv_2   =  zeros(nnzer1(2), 2);
curv_3   =  zeros(nnzer1(3), 2);
curv_4   =  zeros(nnzer1(4), 2);
curv_5   =  zeros(nnzer2(1), 2);
curv_6   =  zeros(nnzer2(2), 2);
curv_7   =  zeros(nnzer2(3), 2);
curv_8   =  zeros(nnzer2(4), 2);
curv_9   =  zeros(nnzer3(1), 2);
curv_10  =  zeros(nnzer3(2), 2);
curv_11  =  zeros(nnzer3(3), 2);
curv_12  =  zeros(nnzer3(4), 2);
curv_13  =  zeros(nnzer4(1), 2);
curv_14  =  zeros(nnzer4(2), 2);
curv_15  =  zeros(nnzer4(3), 2);
curv_16  =  zeros(nnzer4(4), 2);


for k = (Adders_max-nnzer1(1)+1):Adders_max
    curv_1((k+nnzer1(1)-Adders_max), :)  = curva_aux1(k, :, 1);
end
for k = (Adders_max-nnzer1(2)+1):Adders_max
    curv_2((k+nnzer1(2)-Adders_max), :)  = curva_aux1(k, :, 2);
end
for k = (Adders_max-nnzer1(3)+1):Adders_max
    curv_3((k+nnzer1(3)-Adders_max), :)  = curva_aux1(k, :, 3);
end
for k = (Adders_max-nnzer1(4)+1):Adders_max
    curv_4((k+nnzer1(4)-Adders_max), :)  = curva_aux1(k, :, 4);
end
for k = (Adders_max-nnzer2(1)+1):Adders_max
    curv_5((k+nnzer2(1)-Adders_max), :)  = curva_aux2(k, :, 1);
end
for k = (Adders_max-nnzer2(2)+1):Adders_max
    curv_6((k+nnzer2(2)-Adders_max), :)  = curva_aux2(k, :, 2);
end
for k = (Adders_max-nnzer2(3)+1):Adders_max
    curv_7((k+nnzer2(3)-Adders_max), :)  = curva_aux2(k, :, 3);
end
for k = (Adders_max-nnzer2(4)+1):Adders_max
    curv_8((k+nnzer2(4)-Adders_max), :)  = curva_aux2(k, :, 4);
end
for k = (Adders_max-nnzer3(1)+1):Adders_max
    curv_9((k+nnzer3(1)-Adders_max), :)  = curva_aux3(k, :, 1);
end
for k = (Adders_max-nnzer3(2)+1):Adders_max
    curv_10((k+nnzer3(2)-Adders_max), :)  = curva_aux3(k, :, 2);
end
for k = (Adders_max-nnzer3(3)+1):Adders_max
    curv_11((k+nnzer3(3)-Adders_max), :)  = curva_aux3(k, :, 3);
end
for k = (Adders_max-nnzer3(4)+1):Adders_max
    curv_12((k+nnzer3(4)-Adders_max), :)  = curva_aux3(k, :, 4);
end
for k = (Adders_max-nnzer4(1)+1):Adders_max
    curv_13((k+nnzer4(1)-Adders_max), :)  = curva_aux4(k, :, 1);
end
for k = (Adders_max-nnzer4(2)+1):Adders_max
    curv_14((k+nnzer4(2)-Adders_max), :)  = curva_aux4(k, :, 2);
end
for k = (Adders_max-nnzer4(3)+1):Adders_max
    curv_15((k+nnzer4(3)-Adders_max), :)  = curva_aux4(k, :, 3);
end
for k = (Adders_max-nnzer4(4)+1):Adders_max
    curv_16((k+nnzer4(4)-Adders_max), :)  = curva_aux4(k, :, 4);
end





figure(2),...
    p1 = plot(curv_1(:,1),curv_1(:,2)); hold on, ...
    p2 = plot(curv_2(:,1),curv_2(:,2)); ...
    p3 = plot(curv_3(:,1),curv_3(:,2)); ...
    p4 = plot(curv_4(:,1),curv_4(:,2)); ...    
    ax = gca;
    ax.ColorOrderIndex = 1;
    p5 = plot(curv_5(:,1),curv_5(:,2), '--'); ...
    p6 = plot(curv_6(:,1),curv_6(:,2), '--'); ...
    p7 = plot(curv_7(:,1),curv_7(:,2), '--'); ...
    p8 = plot(curv_8(:,1),curv_8(:,2), '--'); grid on, hold off, ...
    legend('Off-line - N = 2', 'Off-line - N = 3', 'Off-line - N = 4', 'Off-line - N = 5', 'In-loop - N = 2', 'In-loop - N = 3', 'In-loop - N = 4', 'In-loop - N = 5'), ...
    xlabel('Adders'), ylabel('MSE [dB]'), ...
    p1.LineWidth = 1;
    p2.LineWidth = 1;
    p3.LineWidth = 1;
    p4.LineWidth = 1;
    p5.LineWidth = 1;
    p6.LineWidth = 1;
    p7.LineWidth = 1;
    p8.LineWidth = 1;


figure(3),...
    p1 = plot(curv_9(:,1),curv_9(:,2)); hold on, ...
    p2 = plot(curv_10(:,1),curv_10(:,2)); ...
    p3 = plot(curv_11(:,1),curv_11(:,2)); ...
    p4 = plot(curv_12(:,1),curv_12(:,2)); ...
    ax = gca;
    ax.ColorOrderIndex = 1;
    p5 = plot(curv_13(:,1),curv_13(:,2), '--'); ...
    p6 = plot(curv_14(:,1),curv_14(:,2), '--'); ...
    p7 = plot(curv_15(:,1),curv_15(:,2), '--'); ...
    p8 = plot(curv_16(:,1),curv_16(:,2), '--'); grid on, hold off, ...
    legend('Off-line - N = 2', 'Off-line - N = 3', 'Off-line - N = 4', 'Off-line - N = 5', 'In-loop - N = 2', 'In-loop - N = 3', 'In-loop - N = 4', 'In-loop - N = 5'), ...
    xlabel('Adders'), ylabel('MSE [dB]'), ...
    p1.LineWidth = 1;
    p2.LineWidth = 1;
    p3.LineWidth = 1;
    p4.LineWidth = 1;
    p5.LineWidth = 1;
    p6.LineWidth = 1;
    p7.LineWidth = 1;
    p8.LineWidth = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Graficos (2) MSE - Adders/Coeff.:

%  MSE - NLMS/SOPOT e NLMS/Inloop-SOPOT

N  = [1 2 3 4];

curva1          = zeros(Adders_max, 2, 4);
curva2          = zeros(Adders_max, 2, 4);
curva3          = zeros(Adders_max, 2, 4);
curva4          = zeros(Adders_max, 2, 4);
curva_aux1      = zeros(Adders_max, 2, 4);
curva_aux2      = zeros(Adders_max, 2, 4);
curva_aux3      = zeros(Adders_max, 2, 4);
curva_aux4      = zeros(Adders_max, 2, 4);
nnzer1          = zeros(1, 4);
nnzer2          = zeros(1, 4);
nnzer3          = zeros(1, 4);
nnzer4          = zeros(1, 4);

tabela_grafic   = sortrows(tabela, 2);

for nAdders = 1:4
    
    curva1(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 3);    
    curva1(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 6);
    
    curva2(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 3);    
    curva2(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 10);
    
    curva3(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 3);    
    curva3(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 7);
    
    curva4(:, 1, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 3);    
    curva4(:, 2, nAdders)     =  tabela_grafic(((Adders_max*(N(nAdders)-1))+1):(Adders_max*N(nAdders)), 8);
    

    for k = 2:Adders_max
        if curva1(k, 2, nAdders) >= curva1(k-1, 2, nAdders)
        else
            curva_aux1(k, :, nAdders) = curva1(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva2(k, 2, nAdders) >= curva2(k-1, 2, nAdders)
        else
            curva_aux2(k, :, nAdders) = curva2(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva3(k, 2, nAdders) >= curva3(k-1, 2, nAdders)
        else
            curva_aux3(k, :, nAdders) = curva3(k, :, nAdders);
        end
    end
    
    for k = 2:Adders_max
        if curva4(k, 2, nAdders) >= curva4(k-1, 2, nAdders)
        else
            curva_aux4(k, :, nAdders) = curva4(k, :, nAdders);
        end
    end   
    
    nnzer1(nAdders)              =   nnz(curva_aux1(:, 1, nAdders));
    nnzer2(nAdders)              =   nnz(curva_aux2(:, 1, nAdders));
    nnzer3(nAdders)              =   nnz(curva_aux3(:, 1, nAdders));
    nnzer4(nAdders)              =   nnz(curva_aux4(:, 1, nAdders));
    curva_aux1(:, :, nAdders)    =   sortrows(curva_aux1(:, :, nAdders), 1);
    curva_aux2(:, :, nAdders)    =   sortrows(curva_aux2(:, :, nAdders), 1);
    curva_aux3(:, :, nAdders)    =   sortrows(curva_aux3(:, :, nAdders), 1);
    curva_aux4(:, :, nAdders)    =   sortrows(curva_aux4(:, :, nAdders), 1);
    

end

curv_1   =  zeros(nnzer1(1), 2);
curv_2   =  zeros(nnzer1(2), 2);
curv_3   =  zeros(nnzer1(3), 2);
curv_4   =  zeros(nnzer1(4), 2);
curv_5   =  zeros(nnzer2(1), 2);
curv_6   =  zeros(nnzer2(2), 2);
curv_7   =  zeros(nnzer2(3), 2);
curv_8   =  zeros(nnzer2(4), 2);
curv_9   =  zeros(nnzer3(1), 2);
curv_10  =  zeros(nnzer3(2), 2);
curv_11  =  zeros(nnzer3(3), 2);
curv_12  =  zeros(nnzer3(4), 2);
curv_13  =  zeros(nnzer4(1), 2);
curv_14  =  zeros(nnzer4(2), 2);
curv_15  =  zeros(nnzer4(3), 2);
curv_16  =  zeros(nnzer4(4), 2);


for k = (Adders_max-nnzer1(1)+1):Adders_max
    curv_1((k+nnzer1(1)-Adders_max), :)  = curva_aux1(k, :, 1);
end
for k = (Adders_max-nnzer1(2)+1):Adders_max
    curv_2((k+nnzer1(2)-Adders_max), :)  = curva_aux1(k, :, 2);
end
for k = (Adders_max-nnzer1(3)+1):Adders_max
    curv_3((k+nnzer1(3)-Adders_max), :)  = curva_aux1(k, :, 3);
end
for k = (Adders_max-nnzer1(4)+1):Adders_max
    curv_4((k+nnzer1(4)-Adders_max), :)  = curva_aux1(k, :, 4);
end
for k = (Adders_max-nnzer2(1)+1):Adders_max
    curv_5((k+nnzer2(1)-Adders_max), :)  = curva_aux2(k, :, 1);
end
for k = (Adders_max-nnzer2(2)+1):Adders_max
    curv_6((k+nnzer2(2)-Adders_max), :)  = curva_aux2(k, :, 2);
end
for k = (Adders_max-nnzer2(3)+1):Adders_max
    curv_7((k+nnzer2(3)-Adders_max), :)  = curva_aux2(k, :, 3);
end
for k = (Adders_max-nnzer2(4)+1):Adders_max
    curv_8((k+nnzer2(4)-Adders_max), :)  = curva_aux2(k, :, 4);
end
for k = (Adders_max-nnzer3(1)+1):Adders_max
    curv_9((k+nnzer3(1)-Adders_max), :)  = curva_aux3(k, :, 1);
end
for k = (Adders_max-nnzer3(2)+1):Adders_max
    curv_10((k+nnzer3(2)-Adders_max), :)  = curva_aux3(k, :, 2);
end
for k = (Adders_max-nnzer3(3)+1):Adders_max
    curv_11((k+nnzer3(3)-Adders_max), :)  = curva_aux3(k, :, 3);
end
for k = (Adders_max-nnzer3(4)+1):Adders_max
    curv_12((k+nnzer3(4)-Adders_max), :)  = curva_aux3(k, :, 4);
end
for k = (Adders_max-nnzer4(1)+1):Adders_max
    curv_13((k+nnzer4(1)-Adders_max), :)  = curva_aux4(k, :, 1);
end
for k = (Adders_max-nnzer4(2)+1):Adders_max
    curv_14((k+nnzer4(2)-Adders_max), :)  = curva_aux4(k, :, 2);
end
for k = (Adders_max-nnzer4(3)+1):Adders_max
    curv_15((k+nnzer4(3)-Adders_max), :)  = curva_aux4(k, :, 3);
end
for k = (Adders_max-nnzer4(4)+1):Adders_max
    curv_16((k+nnzer4(4)-Adders_max), :)  = curva_aux4(k, :, 4);
end



figure(4),
plot(curv_1(:,1),curv_1(:,2)), hold on, ...
    plot(curv_2(:,1),curv_2(:,2)), ...
    plot(curv_3(:,1),curv_3(:,2)), ...
    plot(curv_4(:,1),curv_4(:,2)), ...    
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(curv_5(:,1),curv_5(:,2), '--'), ...
    plot(curv_6(:,1),curv_6(:,2), '--'), ...
    plot(curv_7(:,1),curv_7(:,2), '--'), ...
    plot(curv_8(:,1),curv_8(:,2), '--'), grid on, hold off, ...
    legend('Off-line - N = 2', 'Off-line - N = 3', 'Off-line - N = 4', 'Off-line - N = 5', 'In-loop - N = 2', 'In-loop - N = 3', 'In-loop - N = 4', 'In-loop - N = 5'), ...
    xlabel('Adders/Coeff.'), ylabel('MSE [dB]'), ...


figure(5),...
    plot(curv_9(:,1),curv_9(:,2)), hold on, ...
    plot(curv_10(:,1),curv_10(:,2)), ...
    plot(curv_11(:,1),curv_11(:,2)), ...
    plot(curv_12(:,1),curv_12(:,2)), ...
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(curv_13(:,1),curv_13(:,2), '--'), ...
    plot(curv_14(:,1),curv_14(:,2), '--'), ...
    plot(curv_15(:,1),curv_15(:,2), '--'), ...
    plot(curv_16(:,1),curv_16(:,2), '--'), grid on, hold off, ...
    legend('Off-line - N = 2', 'Offline - N = 3', 'Offline - N = 4', 'Offline - N = 5', 'Inloop - N = 2', 'Inloop - N = 3', 'Inloop - N = 4', 'Inloop - N = 5'), ...
    xlabel('Adders/Coeff.'), ylabel('MSE [dB]'), ...




%   Avaliacao do erro da resposta do inloop e do NLMS
erro = abs(e_spt_aux - e_NLMS_aux);

figure(6),
plot(erro(:,end-3,1)), hold on,...
    plot(erro(:,end-1,1)),...
    plot(erro(:,end,1)), grid on, hold off,...
    legend('Adders = 19', 'Adders = 21', 'Adders = 22'),...
    eixoy = '$|\hat{e}^{2}(k)-e^{2}(k)|$';
    xlabel('Number of Iterations'), ylabel(eixoy , 'Interpreter', 'latex'),...




%   Avaliacao dos coeficientes:
figure(7),
scatter(w_coeff_add_N(1, :, end, 1),w_coeff_add_N(2, :, end, 1), '.'), hold on,...  
    scatter(w_spt_coeff_add_N(1, :, end-16, 1),w_spt_coeff_add_N(2, :, end-16, 1), '.'),...
    scatter(w_spt_coeff_add_N(1, :, end-10, 1),w_spt_coeff_add_N(2, :, end-10, 1), '.'),...   
    scatter(w_spt_coeff_add_N(1, :, end, 1),w_spt_coeff_add_N(2, :, end, 1), '.'),...
    
    scatter(w_coeff_add_N(1, 1, end, 1),w_coeff_add_N(2, 1, end, 1), 80, 'sk', 'filled'),...
    ax = gca;
    ax.ColorOrderIndex = 1;
    s1 = scatter(w_coeff_add_N(1, end, end, 1),w_coeff_add_N(2, end, end, 1), 60, 's', 'filled');... 
    s2 = scatter(w_spt_coeff_add_N(1, end, end-16, 1),w_spt_coeff_add_N(2, end, end-16, 1), 60, 's', 'filled');...
    s3 = scatter(w_spt_coeff_add_N(1, end, end-10, 1),w_spt_coeff_add_N(2, end, end-10, 1), 60, 's', 'filled');...
    s4 = scatter(w_spt_coeff_add_N(1, end, end, 1),w_spt_coeff_add_N(2, end, end, 1), 60, 's', 'filled');...
    
    s1.MarkerEdgeColor = 'k';
    s2.MarkerEdgeColor = 'k';
    s3.MarkerEdgeColor = 'k';
    s4.MarkerEdgeColor = 'k';
  
    xlabel('First Coefficient'), ylabel('Second Coefficient'),...
    grid on,...
    
    xmin = w_coeff_add_N(1, end, end, 1) - 0.05;...
    xmax = w_coeff_add_N(1, end, end, 1) + 0.05;...
    ymin = w_coeff_add_N(2, end, end, 1) - 0.05;...
    ymax = w_coeff_add_N(2, end, end, 1) + 0.05;...
    plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k', 'LineWidth', 0.1),...
    legend('NLMS', 'In-loop 6 Adders', 'In-loop 12 Adders', 'In-loop 22 Adders', 'Starting Point', 'Last NLMS Point', 'Last In-loop 6 Adders Point', 'Last In-loop 12 Adders Point', 'Last In-loop 22 Adders Point', 'Zooming Area'),...
    
    ax2 = axes('position', [0.65 0.65 0.245 0.245]);...
        scatter(w_coeff_add_N(1, :, end, 1),w_coeff_add_N(2, :, end, 1), '.'), hold on,...
        xlim(ax2, [xmin xmax]);...
        ylim(ax2, [ymin ymax]);...
        scatter(w_spt_coeff_add_N(1, :, end-16, 1),w_spt_coeff_add_N(2, :, end-16, 1), '.'),...
        scatter(w_spt_coeff_add_N(1, :, end-10, 1),w_spt_coeff_add_N(2, :, end-10, 1), '.'),...   
        scatter(w_spt_coeff_add_N(1, :, end, 1),w_spt_coeff_add_N(2, :, end, 1), '.'), grid on, grid minor, ...
        ax = gca;
        ax.ColorOrderIndex = 1;
        s1 = scatter(w_coeff_add_N(1, end, end, 1),w_coeff_add_N(2, end, end, 1), 80, 's', 'filled');... 
        s2 = scatter(w_spt_coeff_add_N(1, end, end-16, 1),w_spt_coeff_add_N(2, end, end-16, 1), 70, 's', 'filled');...
        s3 = scatter(w_spt_coeff_add_N(1, end, end-10, 1),w_spt_coeff_add_N(2, end, end-10, 1), 70, 's', 'filled');...
        s4 = scatter(w_spt_coeff_add_N(1, end, end, 1),w_spt_coeff_add_N(2, end, end, 1), 70, 's', 'filled');...
        s1.MarkerEdgeColor = 'k';
        s2.MarkerEdgeColor = 'k';
        s3.MarkerEdgeColor = 'k';
        s4.MarkerEdgeColor = 'k';
        hold off;
    hold off;


figure(8),...
    p1 = plot(w_coeff_add_N(1, :, end, 1)); hold on,...
    p2 = plot(w_spt_coeff_add_N(1, :, end-16, 1));...
    p3 = plot(w_spt_coeff_add_N(1, :, end-10, 1));...
    p4 = plot(w_spt_coeff_add_N(1, :, end, 1));...
    ax = gca;
    ax.ColorOrderIndex = 1;
    p5 = plot(w_coeff_add_N(2, :, end, 1), '--');...
    p6 = plot(w_spt_coeff_add_N(2, :, end-16, 1), '--');...
    p7 = plot(w_spt_coeff_add_N(2, :, end-10, 1), '--');...
    p8 = plot(w_spt_coeff_add_N(2, :, end, 1), '--'); grid on, hold off,...
    legend('NLMS - 1rst Coeff.', '6 Adders - 1rst Coeff.', '12 Adders - 1rst Coeff.', '22 Adders - 1rst Coeff.', 'NLMS - 2nd Coeff.', '6 Adders - 2nd Coeff.', '12 Adders - 2nd Coeff.', '22 Adders - 2nd Coeff.'),...
    xlabel('Number of Iterations'), ylabel('Coefficient Value'),...
    p1.LineWidth = 1;
    p2.LineWidth = 1;
    p3.LineWidth = 1;
    p4.LineWidth = 1;
    p5.LineWidth = 1;
    p6.LineWidth = 1;
    p7.LineWidth = 1;
    p8.LineWidth = 1;


%%%%%%%%%%% SNR:

N = [1 2 3 4];

tabela_SNR = sortrows(tabela_aux3, 2);
curva_SNR  = zeros(Adders_max, 5, N_max-1);

for n = 1:4
    curva_SNR(:, 1, n) = tabela_SNR(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 1);
    curva_SNR(:, 3, n) = tabela_SNR(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 14);
    curva_SNR(:, 2, n) = tabela_SNR(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 12);
    curva_SNR(:, 5, n) = tabela_SNR(((Adders_max*(N(n)-1))+1):(Adders_max*N(n)), 13);
end

figure(9),
    ax = gca;
    p1 = plot(curva_SNR(:,1,1), curva_SNR(:,2,1), '--'); hold on,...
    p2 = plot(curva_SNR(:,1,1), curva_SNR(:,2,2), '--');...
    p3 = plot(curva_SNR(:,1,1), curva_SNR(:,2,3), '--');...
    p4 = plot(curva_SNR(:,1,1), curva_SNR(:,2,4), '--');...
    ax.ColorOrderIndex = 1;
    p5 = plot(curva_SNR(:,1,1), curva_SNR(:,3,1));...
    p6 = plot(curva_SNR(:,1,1), curva_SNR(:,3,2));...
    p7 = plot(curva_SNR(:,1,1), curva_SNR(:,3,3));...
    p8 = plot(curva_SNR(:,1,1), curva_SNR(:,3,4));...
    ax.ColorOrderIndex = 1;
    p9 = plot(curva_SNR(:,1,1), curva_SNR(:,5,1), '-.');...
    p10 = plot(curva_SNR(:,1,1), curva_SNR(:,5,2), '-.');...
    p11 = plot(curva_SNR(:,1,1), curva_SNR(:,5,3), '-.');...
    p12 = plot(curva_SNR(:,1,1), curva_SNR(:,5,4), '-.');grid on, hold off,...
    legend('NLMS - N = 2', 'NLMS - N = 3', 'NLMS - N = 4', 'NLMS - N = 5', 'In-loop - N = 2', 'In-loop - N = 3', 'In-loop - N = 4', 'In-loop - N = 5', 'Off-line - N = 2', 'Off-line - N = 3', 'Off-line - N = 4', 'Off-line - N = 5'),...
    xlabel('Adders'), ylabel('SNR [dB]'),...
    p1.LineWidth = 1;
    p2.LineWidth = 1;
    p3.LineWidth = 1;
    p4.LineWidth = 1;
    p5.LineWidth = 1;
    p6.LineWidth = 1;
    p7.LineWidth = 1;
    p8.LineWidth = 1;
    p9.LineWidth = 1;
    p10.LineWidth = 1;
    p11.LineWidth = 1;
    p12.LineWidth = 1;


% Montar comparativo entre o efeito de anulacao da AR de 1a ordem do NLMS e
% da aproximacao por SPT freqz(wn, a)