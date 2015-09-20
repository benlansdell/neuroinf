%Set wd to working directory
wd = './results4_gauss_move_5hz_maxit20_5frame';
%wd = '/Users/dk/Dropbox/AFK_Neuron/Code_GLM_ben/monkeyresults/coherence_networkglm/';
load([wd '/preprocessed_networkglm_sims.mat'])

%For each unit
for icell = 1:nU
    truesp = {};
    simsp = {};
    simsp_cpl = {};
    Coh_cpl = [];
    Coh = [];
    Power = [];
    Power_cpl = [];
    Power_true = [];
    runs = size(proc_withheld.trialstartend,1);
    K = [];
    N_sample_max = 0;
    uncoupled = load([wd '/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
    nspks = 0;
    for i = 1:runs
        ts = proc_withheld.trialstartend(i, 1);
        te = proc_withheld.trialstartend(i, 2);
        truesp{i} = proc_withheld.spiketrain(:,goodunits(icell));
        truesp{i} = truesp{i}(ts:te);
        simsp_cpl{i} = Rt_glm{icell}(:);
        simsp_cpl{i} = simsp_cpl{i}(ts:te);
        simsp{i} = uncoupled.Rt_glm(:);
        simsp{i} = simsp{i}(ts:te);
        N_sample_max = max(te-ts+1, N_sample_max);
        nspks = nspks + sum(truesp{i});
    end
    Pad = 2^(1+nextpow2(N_sample_max));  
    f=linspace(0, 1, Pad)*RefreshRate;

    for i = 1:runs
        %Compute coherence of all runs
        [coh_cpl, pwr_true, pwr_cpl, k] = coherence(truesp{i}, simsp_cpl{i}, N_sample_max);
        [coh, pwr1, pwr]  = coherence(truesp{i}, simsp{i}, N_sample_max);
        Coh_cpl(i,:) = coh_cpl;
        Coh(i,:) = coh;
        Power(i,:) = pwr;
        Power_cpl(i,:) = pwr_cpl;
        Power_true(i,:) = pwr_true;
        K(i) = k;
    end

    %Only compute for not NaN results
    Coh = Coh(~isnan(Coh(:,1)),:);
    size(Coh);
    if size(Coh,1) > 1
        Coh = mean(Coh);
    end
    Coh_cpl = Coh_cpl(~isnan(Coh_cpl(:,1)),:);
    size(Coh_cpl);
    if size(Coh_cpl,1) > 1
        Coh_cpl = mean(Coh_cpl);
    end
    K = mean(K(~isnan(Coh(:,1))));
    DF = runs*K-1; Null = sqrt(1-0.05^(1/DF));
    L = 5; 
    ii = 1:fix(Pad/L);
    subplot(10,1,icell)
    [ax, h1, h2] = plotyy(f(ii), [abs(Coh(ii)); abs(Coh_cpl(ii)); Null*ones(size(ii))], f(ii), [Power(1,ii); Power_cpl(1,ii)]);
    set(h2, 'LineStyle', '--');
    legend('Uncoupled', 'Coupled')
    title(['Unit: ' num2str(goodunits(icell)) ' no. spikes: ' num2str(nspks)])
    ylabel('<|coh|>')
end
saveplot(gcf, [wd '/GLM_coherence_' num2str(goodunits(icell)) '_inmovement.eps'])
