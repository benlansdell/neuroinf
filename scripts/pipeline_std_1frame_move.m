%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nU = 45;                        %no. units
nS = 4;                         %no. stim components
frames = 1;                     %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 20;                      %no. sim repetitions
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
nB = size(proc.stim, 1);
fn_out = './results_std_1frame_move/';
trim = 1;
pca = 0;
dt_glm = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs = {};
for icell = 1:nU
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
    gg0 = makeFittingStruct_GLM_monkey(sta,dt_glm);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', 10};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
    ggs{icell} = gg;
    gg.ih(end-1) = 100;
    %After fitting, put abs refractory period down lower...

    %Simulation with test stim
    disp(num2str(icell));
    Tt = size(proc_withheld.stim,1);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, proc_withheld.stim/p, time_limit);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + conv;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_cell_' num2str(icell) '.mat'],...
        'gg'); %, 'Rt_glm');
end

%Save all
save([fn_out '/all_units.mat'], 'ggs');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 Plot uncoupled filters%
%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filters(ggs, proc, [fn_out '/all_units_filters.eps']);

%%%%%%%%%%%%%%%%%
%Simulate models%
%%%%%%%%%%%%%%%%%

time_limit = 40;
units_conv = zeros(nU,1);
for icell = 1:nU
    load([fn_out '/GLM_cell_' num2str(icell) '.mat']);
    %Simulation with test stim
    disp(num2str(icell));
    stim = proc_withheld.stim;
    stim = stim/p;
    Tt = size(proc_withheld.stim,1);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, stim, time_limit);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + conv;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot uncoupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icell = 1:nU
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
    plot(tidx*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
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
    plot(tidx*proc.binsize, 500*proc_withheld.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc_withheld.cursor(tidx,3), 'r');
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

goodunits = [3 4 6 7 8 9 14 15 16 17 18 19 21 23 27 31 33 34 35 36 37 39 40 41];
nU = length(goodunits);
%Remove the 'bad' units from the dataset
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

%Run fitting...
ggs_cpl = {};
%For testing
maxiter = 10;
%For reals
maxiter = 30;
dt_glm = 0.01;
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
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec(sta,dt_glm,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);
    gg0.ihbas2 = gg0.ihbas;

    opts = {'display', 'iter', 'maxiter', maxiter};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
    ggs_cpl{icell} = gg;

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
    save([fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'],...
        'gg'); %, 'Rt_glm');
end

%Save all
save([fn_out '/all_units_network.mat'], 'ggs_cpl');

%Plot
plot_filters_network(ggs_cpl, proc, [fn_out '/all_units_network_filters.eps']);

load([fn_out '/all_units_network.mat'])
time_limit = 30;
nU = length(ggs_cpl);
%Simulate network model...
stim = proc_withheld.stim;
stim = stim/p;
clear proc;
clear proc_withheld;
%Make sim struct
glm_cpl = makeSimStruct_GLMcpl(ggs_cpl{:});
%Simulation with test stim
Tt = size(stim,1);
Rt_glm = {};
for idx = 1:nU
    Rt_glm{idx} = zeros(Tt,1);
end

for ir = 1:nRep
    ir
    for idx = 1:nU
        Rt_glm{idx}(ceil(iR_glm{idx})) = Rt_glm{idx}(ceil(iR_glm{idx}))+1;
    end
end
for idx = 1:nU
    Rt_glm{idx} = Rt_glm{idx}'/nRep + 1e-8;
end
save([fn_out '/GLM_coupled_simulation_all.mat'], 'Rt_glm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot coupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fn_out '/all_units_network.mat'])
load([fn_out '/GLM_coupled_simulation_all.mat'])
goodunits = [3 4 6 7 8 9 14 15 16 17 18 19 21 23 27 31 33 34 35 36 37 39 40 41];
nU = length(goodunits);
%Remove the 'bad' units from the dataset
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

for icell = 1:nU
    clf
    sigma_fr = .25;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    tstart = ceil(600*RefreshRate);
    tend = ceil(750*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm{icell}(tidx);
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
    plot(tidx*proc.binsize, gftruesp, tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_' num2str(goodunits(icell)) '_sim.eps'], 'eps', [6 6]);
    
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
    pillowsimsp = Rt_glm{icell}(tidx);
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
    spidx = truesp==1;
    plot(tidx(spidx)*proc.binsize, truesp(spidx)-.95, '.', tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_coupled_' num2str(goodunits(icell)) '_sim_zoom.eps'], 'eps', [6 6]);
end