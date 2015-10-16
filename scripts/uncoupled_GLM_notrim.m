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
nU = 45;                        %no. units
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 747;                      %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);    
nB = size(proc.stim, 1);
fn_out = '/results_notrim/';
trim = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
offset = 1;
mkdir([wd fn_out]);

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs = {};
for icell = goodunits
    disp(['Fitting unit ' num2str(icell)]);
    stim = proc.stim;
    stim = stim/p;
    resp = proc.spikes{icell};
    sptrain = proc.spiketrain(:,icell);

    stacked = proc.stacked;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    sta = reshape(sta,nF,[]);

    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey(sta,dt_glm,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', maxit};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, offset);
    ggs{icell} = gg;

    save([wd fn_out '/GLM_cell_' num2str(icell) '.mat'], 'gg');
end

%Save all
save([wd fn_out '/all_units.mat'], 'ggs');
load([wd fn_out '/all_units.mat']);

%%%%%%%%%%%%%%%%%%%
%3 Simulate models%
%%%%%%%%%%%%%%%%%%%

time_limit = 80;
units_conv = zeros(nU,1);
logl_glm_uncoupled = [];
rng('shuffle')
for icell = goodunits
    load([wd fn_out '/GLM_cell_' num2str(icell) '.mat']);
    %Simulation with test stim
    disp(num2str(icell));
    %Only within trial times...
    Rt = proc_withheld.spiketrain(:,icell);
    stim = proc_withheld.stim(:,:);
    stim = stim/p;
    stim = stim(:,:);
    Tt = size(stim,1);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;

    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, converged] = simGLM_monkey(gg, stim, time_limit, 1);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + converged;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    logl_glm = mean(Rt.*log(Rt_glm)-(Rt_glm)/RefreshRate) ;
    logl_glm_uncoupled(icell) = logl_glm;
    save([wd fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged', 'logl_glm');
end

%%%%%%%%%
%4 Plot %
%%%%%%%%%
plot_filters(ggs, proc, [wd fn_out '/all_units_filters.eps'], goodunits);

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
    load([wd fn_out '/GLM_cell_' num2str(icell) '.mat'])
    load([wd fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'])
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
    saveplot(gcf, [wd fn_out '/GLM_cell_' num2str(icell) '_sim.eps'], 'eps', [6 6]);
end