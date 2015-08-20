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
frames = 100;                   %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 20;                      %no. sim repetitions
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
nB = size(proc.stim, 1);
fn_out = './results_gauss_move/';
trim = 1;
pca = 0;
Dt = 10;

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
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec(sta,dt,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', 100};
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
plot_filters(ggs, proc, [fn_out '/all_units_filters.eps']);

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
    tend = ceil(900*RefreshRate);
    fn_in = [fn_out '/GLM_cell_' num2str(icell) '.mat'];
    load(fn_in,'gg','Rt_glm');
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
    truesp = proc_withheld.spikes(tidx,icell);
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
maxiter = 3;
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
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec(sta,dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);

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