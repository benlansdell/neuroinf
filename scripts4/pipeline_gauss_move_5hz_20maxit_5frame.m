%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
goodunits = [4,7,8,14,15,16,17,20,24,36,37,41];

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nU = 45;                        %no. units
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 20;                      %no. sim repetitions
std = 0;
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);    
nB = size(proc.stim, 1);
fn_out = './results4_gauss_move_5hz_maxit20_5frame/';
trim = 1;
pca = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs = {};
for icell = goodunits
    disp(num2str(icell));
    stim = proc.stim;
    stim = stim/p;
    resp = proc.spikes{icell};
    sptrain = proc.spiketrain(:,icell);

    stacked = proc.stacked;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    sta = reshape(sta,nF,[]);

    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec(sta,dt_glm,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', maxit};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
    ggs{icell} = gg;

    %Simulation with test stim
    %disp(num2str(icell));
    %stim = proc_withheld.stim;
    %stim = stim/p;
    %Tt = size(proc_withheld.stim,1);
    %Rt_glm = zeros(1,Tt);
    %for ir = 1:nRep
    %    ir
    %    [iR_glm,vmem,Ispk] = simGLM_monkey(ggs{icell}, stim);
    %    Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    %end
    %Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_cell_' num2str(icell) '.mat'],...
        'gg'); %, 'Rt_glm');
end

%Save all
save([fn_out '/all_units.mat'], 'ggs');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 Plot uncoupled filters%
%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filters(ggs, proc, [fn_out '/all_units_filters.eps'], goodunits);

%%%%%%%%%%%%%%%%%
%Simulate models%
%%%%%%%%%%%%%%%%%

time_limit = 40;
units_conv = zeros(nU,1);
for icell = goodunits
    load([fn_out '/GLM_cell_' num2str(icell) '.mat']);
    %Simulation with test stim
    disp(num2str(icell));
    Tt = size(proc_withheld.stim,1);
    Rt = proc_withheld.spiketrain(:,icell);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, proc_withheld.stim/p, time_limit, 0);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + conv;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    logl_glm = mean(Rt.*log(Rt_glm)-(Rt_glm)*dt) ;
    save([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged', 'logl_glm');
end

logl_glm_uncoupled = [];
for i = goodunits
    load([fn_out '/GLM_cell_simulation_' num2str(i) '.mat']);
    logl_glm_uncoupled(i) = logl_glm;
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot uncoupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icell = goodunits
    clf
    sigma_fr = .25;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(600*RefreshRate);
    tend = ceil(750*RefreshRate);
    fn_in = [fn_out '/GLM_cell_' num2str(icell) '.mat'];
    load([fn_out '/GLM_cell_' num2str(icell) '.mat'])
    load([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'])
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    plot(tidx*proc.binsize, gftruesp, tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_cell_' num2str(icell) '_sim.eps'], 'eps', [6 6]);
    
    clf
    sigma_fr = .01;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(600*RefreshRate);
    tend = ceil(605*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    spidx = truesp==1;
    plot(tidx(spidx)*proc.binsize, truesp(spidx)-.95, '.', tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_cell_' num2str(icell) '_sim_zoom.eps'], 'eps', [6 6]);
end


%%%%%%%%%%%%%%%%%%%
%Fit network model%
%%%%%%%%%%%%%%%%%%%

goodunits = [4,7,8,14,15,17,20,24,36,41];
nU = length(goodunits);
%Remove the 'bad' units from the dataset
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

%Run fitting...
ggs_cpl = {};
%For testing
maxiter = 10;
%For reals
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
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
    ggs_cpl{icell} = gg;
    save([fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'],...
        'gg'); %, 'Rt_glm');
end

for icell = 1:10
    ggs_cpl{icell}.ihbas2 = ggs_cpl{icell}.ihbas;
end

%Save all
save([fn_out '/all_units_network.mat'], 'ggs_cpl');

%Plot
plot_filters_network_all(ggs_cpl, proc, [fn_out '/all_units_network_filters.eps'], goodunits);
plot_filters_network_compare(ggs_cpl, ggs, proc, [fn_out '/all_units_network_filters_compare.eps'], goodunits);

load([fn_out '/all_units_network.mat'])
nU = length(ggs_cpl);
%Simulate network model...
stim = proc_withheld.stim;
stim = stim/p;
time_limit = 2400;
nRep = 10;
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std); 
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);   
clear proc;
clear proc_withheld;
for i = 1:nU
    ggs_cpl{i}.ihbas2 = ggs_cpl{i}.ihbas;
end
simstruct = makeSimStruct_GLMcpl(ggs_cpl{:});
%Simulation with test stim
disp(num2str(icell));
Tt = size(stim,1);
Rt_glm = {};
for i = 1:nU
    Rt_glm{i} = zeros(1,Tt);
end
for ir = 1:nRep
    ir
    [iR_glm,vmem,Ispk] = simGLM_monkey(simstruct, stim, time_limit);
    for i = 1:nU
        Rt_glm{i}(ceil(iR_glm{i})) = Rt_glm{i}(ceil(iR_glm{i}))+1;
    end
end

logl_glm = [];
for i = 1:nU
    Rt = proc_withheld.spiketrain(:,i);
    Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
    if size(Rt_glm{i},1)==1
        Rt_glm{i} = Rt_glm{i}';
    end
    %Compute log-likelihood:
    logl_glm(i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(dt)) ;
end
save([fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');

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
    tend = ceil(650*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm{idx}(tidx);
    uncoupled = load([fn_out '/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations, chopping out no trial times%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fn_out '/GLM_coupled_simulation.mat'])
for idx = 1:nU
    icell = (idx);
    clf
    sigma_fr = .1;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(200*RefreshRate);
    tend = ceil(550*RefreshRate);
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
    saveplot(gcf, [fn_out '/GLM_coupled_cell_' num2str(goodunits(icell)) '_sim_inmovement.eps'], 'eps', [6 6]);
    
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
    save([fn_out '/GLM_sims_cell_' num2str(goodunits(icell)) '_inmovement.mat'], 'truesp', 'simsp', 'simsp_cpl');
end

%Trial start to trial end
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
