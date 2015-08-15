%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
dt = 0.1;                       %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nU = 45;                        %no. units
nS = 4;                         %no. stim components
nRep = 20;                      %no. sim repetitions
nC = 10;                        %no. stim PC components to keep
p = nC*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, nS);    

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
    sta = stacked'*sptrain/sum(sptrain) - mean(stacked,1)'; 
    sta = reshape(sta,frames,[]);

    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey(sta,dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', 100};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc);

    %Simulation with test stim
    stim = proc_withheld.stim;
    stim = stim/p;
    Tt = size(proc_withheld.stim,1);
    Rt_glm = zeros(1,Tt);
    for ir = 1:nRep
        [iR_glm,vmem,Ispk] = simGLM(gg, stim);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    end
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save(['./results/Uncoupled_GLM_cell_' num2str(icell) '.mat'],...
        'gg', 'Rt_glm');
    ggs{icell} = gg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 Plot uncoupled filters%
%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filters(ggs, proc, fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot uncoupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
sigma_fr = .25;
sigma_fr = sigma_fr*RefreshRate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

tstart = 8001;
tend = 12000;
unitidx = 1;
unit = t(unitidx,2);
load(['GLM_model_' num2str(unit) '_ext_5Hz_trim.mat'],'gg','Rt_glm');
tidx = tstart:tend;
truesp = proc_withheld.spikes(tidx,unitidx);
pillowsimsp = Rt_glm(tidx);
mysimsp = simtrains{unitidx}(tidx);
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
gfmysimsp = conv(mysimsp, gaussFilter_fr, 'same');
plot(tidx*proc.binsize, gftruesp, tidx*proc.binsize, gfpillowsimsp, tidx*proc.binsize, gfmysimsp);
title('n rep: 20')
xlabel('seconds')
ylabel('predicted probability spiking')
legend('true spike train', 'Pillow''s GLM', 'Indep GLM')
saveplot(gcf, [fn_out '_sim'], 'eps', [6 6]);

clf
sigma_fr = .01;
sigma_fr = sigma_fr*RefreshRate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

tstart = 800;
tend = 1200;
unitidx = 7;
unit = t(unitidx,2);
load(['GLM_model_' num2str(unit) '_ext_5Hz_trim.mat'],'gg','Rt_glm');
tidx = tstart:tend;
truesp = proc_withheld.spikes(tidx,unitidx);
pillowsimsp = Rt_glm(tidx);
mysimsp = simtrains{unitidx}(tidx);
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
gfmysimsp = conv(mysimsp, gaussFilter_fr, 'same');
spidx = truesp==1;
plot(tidx(spidx)*proc.binsize, truesp(spidx)-.95, '.', tidx*proc.binsize, gfpillowsimsp, tidx*proc.binsize, gfmysimsp);
title('n rep: 20')
xlabel('seconds')
ylabel('predicted probability spiking')
legend('true spike train', 'Pillow''s GLM', 'Indep GLM')
saveplot(gcf, [fn_out '_sim_zoom'], 'eps', [6 6]);