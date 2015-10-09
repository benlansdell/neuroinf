%Set to working directory
wd = '.';

%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
goodunits = [4,7,14,15,17,20,24,36,41];

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
nRep = 83;                      %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    
nB = size(proc.stim, 1);
fn_out = '/results/';
trim = 1;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir([wd fn_out]);

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting coupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs_cpl = {};
maxiter = 20;
for icell = 1:nU
    disp(['Fitting unit ' num2str(goodunits(icell))]);
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
    gg0 = makeFittingStruct_GLM_monkey(sta,dt,Dt);
    gg0.ihbas2 = gg0.ihbas;
    gg0.tsp = resp';
    gg0.tspi = 1;
    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);
    opts = {'display', 'iter', 'maxiter', maxiter};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim);
    ggs_cpl{icell} = gg;
    save([wd fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'], 'gg');
end

%Save all
save([wd fn_out '/all_units_network.mat'], 'ggs_cpl');
load([wd fn_out '/all_units_network.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate network model%
%%%%%%%%%%%%%%%%%%%%%%%%

stim = proc_withheld.stim;
stim = stim/p;
stim = stim(:,:);
time_limit = 2400;
Tt = size(stim,1);
Rt_glm = {};
for i = 1:nU
    Rt_glm{i} = zeros(1,Tt);
    ggs_cpl{i}.ihbas2 = ggs_cpl{i}.ihbas;
end
rng('shuffle')
simstruct = makeSimStruct_GLMcpl(ggs_cpl{:});
for ir = 1:nRep
    ir
    [iR_glm,vmem,Ispk] = simGLM_monkey(simstruct, stim, time_limit);
    for i = 1:nU
        Rt_glm{i}(ceil(iR_glm{i})) = Rt_glm{i}(ceil(iR_glm{i}))+1;
    end
end

save([wd fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm');

%Compute likelihood
logl_glm = [];
for i = 1:nU
    Rt = proc_withheld.spiketrain(:,i);
    Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
    if size(Rt_glm{i},1)==1
        Rt_glm{i} = Rt_glm{i}';
    end
    %Compute log-likelihood:
    logl_glm(i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
end
save([wd fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');

%Save stuff needed to run this chunk of code
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);    
save([wd fn_out '/preprocessed_networkglm_sims.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')

%Compare likelihood to uncoupled likelihood:
load([wd fn_out '/GLM_coupled_simulation.mat'])
logl_glm_uncoupled = [];
for idx = 1:nU
    icell = goodunits(idx);
    uncoupled = load([wd fn_out '/GLM_cell_simulation_' num2str(icell) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end
clf
scatter(logl_glm_uncoupled, logl_glm)
hold on; plot([-1.5 0], [-1.5 0], 'r')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare.eps'])

%%%%%%%%%%%%%%%%%%%%%%
%Plot network filters%
%%%%%%%%%%%%%%%%%%%%%%

plot_filters_network_all(ggs_cpl, proc, [wd fn_out '/all_units_network_filters.eps'], goodunits);
load([wd fn_out '/all_units.mat'])
plot_filters_network_compare(ggs_cpl, ggs, proc, [wd fn_out '/all_units_network_filters_compare.eps'], goodunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([wd fn_out '/GLM_coupled_simulation.mat'])
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
    uncoupled = load([wd fn_out '/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
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
    saveplot(gcf, [wd fn_out '/GLM_coupled_cell_' num2str(goodunits(icell)) '_sim.eps'], 'eps', [6 6]);

    simsp = uncoupled.Rt_glm;
    simsp_cpl = Rt_glm{idx};
    truesp = proc_withheld.spiketrain(:,icell);
    save([wd fn_out '/GLM_sims_cell_' num2str(goodunits(icell)) '.mat'], 'truesp', 'simsp', 'simsp_cpl');
end
