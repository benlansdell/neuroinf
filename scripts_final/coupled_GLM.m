%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
goodunits = [4,7,8,14,15,17,20,24,36,41];

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nU = length(goodunits);
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 10;                      %no. sim repetitions
std = 0;
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);    
nB = size(proc.stim, 1);
fn_out = './results_coupled_glm/';
trim = 1;
pca = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir(fn_out);

%Remove the 'bad' units from the dataset
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

%Run fitting...
ggs_cpl = {};
maxiter = 20;
for icell = 1:nU
    disp(num2str(icell));
    stim = proc.stim;
    stim = stim/p;
    resp = proc.spikes{icell};
    sptrain = proc.spiketrain(:,icell);
    nicell = [(1:icell-1), (icell+1:nU)];
    %Add coupling to the other spike trains
    coupled = proc.spikes(nicell);
    for idx = 1:length(coupled)
        coupled{idx} = coupled{idx}';
    end
    stacked = proc.stacked;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    sta = reshape(sta,nF,[]);
    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec_refract(sta,dt,Dt);
    gg0.ihbas2 = gg0.ihbas;
    gg0.tsp = resp';
    gg0.tspi = 1;
    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);
    opts = {'display', 'iter', 'maxiter', maxiter};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim,pca);
    ggs_cpl{icell} = gg;
    save([fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'], 'gg');
end

%Save all
save([fn_out '/all_units_network.mat'], 'ggs_cpl');
load([fn_out '/all_units_network.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate network model%
%%%%%%%%%%%%%%%%%%%%%%%%

stim = proc_withheld.stim;
stim = stim/p;
stim = stim(1:20000,:);
time_limit = 2400;
clear proc;
clear proc_withheld;
%Simulation with test stim
Tt = size(stim,1);
Rt_glm = {};
for i = 1:nU
    Rt_glm{i} = zeros(1,Tt);
    ggs_cpl{i}.ihbas2 = ggs_cpl{i}.ihbas;
end

simstruct = makeSimStruct_GLMcpl(ggs_cpl{:});
for ir = 1:nRep
    ir
    [iR_glm,vmem,Ispk] = simGLM_monkey(simstruct, stim, time_limit);
    for i = 1:nU
        Rt_glm{i}(ceil(iR_glm{i})) = Rt_glm{i}(ceil(iR_glm{i}))+1;
    end
end

nRep = 10;
logl_glm = [];
for i = 1:nU
    Rt = proc_withheld.spiketrain(1:20000,i);
    Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
    if size(Rt_glm{i},1)==1
        Rt_glm{i} = Rt_glm{i}';
    end
    %Compute log-likelihood:
    logl_glm(i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
end
save([fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');

%Save stuff needed to run this chunk of code
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);    
save([fn_out '/preprocessed_networkglm_sims.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')

logl_glm_uncoupled = [];
for idx = 1:nU
    icell = goodunits(idx);
    uncoupled = load(['./results_uncoupled_GLM/GLM_cell_simulation_' num2str(icell) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end

%Compare likelihood to uncoupled likelihood:
scatter(logl_glm, logl_glm_uncoupled)
hold on; plot([-1.5 0], [-1.5 0], 'r')
xlabel('Coupled log-likelihood');
ylabel('Uncoupled log-likelihood')
saveplot(gcf, [fn_out '/GLM_loglikelihood_compare.eps'])

%%%%%%%%%%%%%%%%%%%%%%
%Plot network filters%
%%%%%%%%%%%%%%%%%%%%%%

plot_filters_network_all(ggs_cpl, proc, [fn_out '/all_units_network_filters.eps'], goodunits);
load('./results_uncoupled_GLM/all_units.mat')
plot_filters_network_compare(ggs_cpl, ggs, proc, [fn_out '/all_units_network_filters_compare.eps'], goodunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fn_out '/GLM_coupled_simulation.mat'])
for idx = 1:nU
    icell = (idx);
    clf
    sigma_fr = .25;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(500*RefreshRate);
    tend = ceil(550*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm{idx}(tidx);
    uncoupled = load(['./results_uncoupled_GLM/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
    pillowsimsp_uncoupled = uncoupled.Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    gfpillowsimsp_uncoupled = conv(pillowsimsp_uncoupled, gaussFilter_fr, 'same');
    plot(tidx*proc.binsize, gftruesp, tidx*proc.binsize, gfpillowsimsp, tidx*proc.binsize, gfpillowsimsp_uncoupled);
    title(['no. rep: 20. unit: ' num2str(goodunits(idx)) '(' proc.unitnames{idx} ')'])
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Coupled GLM', 'Uncoupled GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_cell_' num2str(goodunits(icell)) '_sim.eps'], 'eps', [6 6]);
    
    clf
    sigma_fr = .01;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(500*RefreshRate);
    tend = ceil(505*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm{idx}(tidx);
    pillowsimsp_uncoupled = uncoupled.Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    gfpillowsimsp_uncoupled = conv(pillowsimsp_uncoupled, gaussFilter_fr, 'same');
    spidx = truesp==1;
    plot(tidx(spidx)*proc.binsize, truesp(spidx)-.95, '.', tidx*proc.binsize, gfpillowsimsp, tidx*proc.binsize, gfpillowsimsp_uncoupled);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Coupled GLM', 'Uncoupled GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_cell_' num2str(goodunits(icell)) '_sim_zoom.eps'], 'eps', [6 6]);

    simsp = uncoupled.Rt_glm;
    simsp_cpl = Rt_glm{idx};
    truesp = proc_withheld.spiketrain(:,icell);
    save([fn_out '/GLM_sims_cell_' num2str(goodunits(icell)) '.mat'], 'truesp', 'simsp', 'simsp_cpl');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations, chopping out no movement times%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fn_out '/GLM_coupled_simulation.mat'])
for idx = 1:nU
    clf
    sigma_fr = .1;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(1*RefreshRate);
    tend = ceil(200*RefreshRate);
    nB = size(proc_withheld.stim,1);
    tidx = ((1:nB)>tstart) & ((1:nB)<tend) & (proc_withheld.intrial==1)';
    tt = 1:sum(tidx == 1);
    uncoupled = load(['./results_uncoupled_GLM/GLM_cell_simulation_' num2str(goodunits(idx)) '.mat']);
    trialstart = [0; (diff(proc_withheld.intrial) == 1)];
    gftruesp = conv(proc_withheld.spiketrain(:,idx), gaussFilter_fr, 'same');
    gfpillowsimsp = conv(Rt_glm{idx}(:), gaussFilter_fr, 'same');
    gfpillowsimsp_uncoupled = conv(uncoupled.Rt_glm(:), gaussFilter_fr, 'same');
    gftruesp = gftruesp(tidx);
    gfpillowsimsp = gfpillowsimsp(tidx);
    gfpillowsimsp_uncoupled = gfpillowsimsp_uncoupled(tidx);
    trialstart = trialstart(tidx);
    trialstart = find(trialstart);

    subplot(2,1,1)
    hold on
    plot(tt*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
    hold on
    for ts = trialstart
        plot([ts, ts]*proc.binsize, [0 1], 'Color', [0.3 0.3 0.3])
    end
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    plot(tt*proc.binsize, gftruesp, tt*proc.binsize, gfpillowsimsp, tt*proc.binsize, gfpillowsimsp_uncoupled);
    hold on
    for ts = trialstart
        plot([ts, ts]*proc.binsize, [0 1], 'Color', [0.3 0.3 0.3])
    end
    title(['no. rep: 20. unit: ' num2str(goodunits(idx)) '(' proc.unitnames{idx} ')'])
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Coupled GLM', 'Uncoupled GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_cell_' num2str(goodunits(idx)) '_sim_inmovement.eps'], 'eps', [6 6]);
    save([fn_out '/GLM_sims_cell_' num2str(goodunits(idx)) '_inmovement.mat'], 'truesp', 'simsp', 'simsp_cpl');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations, chopping out no trial times%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames);  
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

for idx = 1:nU
    icell = (idx);
    clf
    sigma_fr = .1;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(500*RefreshRate);
    tend = ceil(650*RefreshRate);
    nB = size(proc_withheld.stim,1);
    tidx = ((1:nB)>tstart) & ((1:nB)<tend) & (proc_withheld.intrial==1)';
    tt = 1:sum(tidx == 1);
    uncoupled = load([fn_out '/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
    trialstart = [0; (diff(proc_withheld.intrial) == 1)];
    gftruesp = conv(proc_withheld.spiketrain(:,icell), gaussFilter_fr, 'same');
    gfpillowsimsp = conv(Rt_glm{idx}(:), gaussFilter_fr, 'same');
    gfpillowsimsp_uncoupled = conv(uncoupled.Rt_glm(:), gaussFilter_fr, 'same');
    gftruesp = gftruesp(tidx);
    gfpillowsimsp = gfpillowsimsp(tidx);
    gfpillowsimsp_uncoupled = gfpillowsimsp_uncoupled(tidx);
    trialstart = trialstart(tidx);
    trialstart = find(trialstart);

    subplot(2,1,1)
    hold on
    plot(tt*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tt*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
    hold on
    for ts = trialstart
        plot([ts, ts]*proc.binsize, [0 1], 'Color', [0.3 0.3 0.3])
    end
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    plot(tt*proc.binsize, gftruesp, tt*proc.binsize, gfpillowsimsp, tt*proc.binsize, gfpillowsimsp_uncoupled);
    hold on
    for ts = trialstart
        plot([ts, ts]*proc.binsize, [0 1], 'Color', [0.3 0.3 0.3])
    end
    title(['no. rep: 20. unit: ' num2str(goodunits(idx)) '(' proc.unitnames{idx} ')'])
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Coupled GLM', 'Uncoupled GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_cell_' num2str(goodunits(icell)) '_sim_intrial.eps'], 'eps', [6 6]);

    truesp = {};
    simsp = {};
    simsp_cpl = {};
    for i = 1:size(proc_withheld.trialstartend,1)
        ts = proc_withheld.trialstartend(i, 1);
        te = proc_withheld.trialstartend(i, 2);
        truesp{i} = proc_withheld.spiketrain(:,icell);
        truesp{i} = truesp{i}(ts:te);
        simsp_cpl{i} = Rt_glm{icell}(:);
        simsp_cpl{i} = simsp_cpl{i}(ts:te);
        simsp{i} = uncoupled.Rt_glm(:);
        simsp{i} = simsp{i}(ts:te);
    end
    save([fn_out '/GLM_sims_cell_' num2str(goodunits(icell)) '_intrial.mat'], 'truesp', 'simsp', 'simsp_cpl');
end