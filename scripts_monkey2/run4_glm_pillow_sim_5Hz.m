global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 5; 
dt = 0.005;
datafile = './data/mabel_reaching_5-4-10.mat';
%Test run...
%nU = 4;
nU = 45;
nS = 4;
samplerate = RefreshRate;
nRep = 20;
frames = 9;
binsize = 1/RefreshRate;
processed = preprocess_monkey_pillow_5Hz(datafile, binsize, dt, frames);    

sigma_fr = .25;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

nB = size(processed.spiketrain,1);
trainmax = nB*0.8;
%Split into 20% test and 80% training set
processedtest = processed;
trainingidx = (1:nB)<=trainmax;
processed.cursor = processed.cursor(trainingidx,:);
processed.grip = processed.grip(trainingidx,:);
processed.stim = processed.stim(trainingidx,:);
processed.spiketrain = processed.spiketrain(trainingidx,:);
processed.stacked = processed.stacked(trainingidx,:);
processed.intrial = processed.intrial(trainingidx);

processedtest.cursor = processedtest.cursor(~trainingidx,:);
processedtest.grip = processedtest.grip(~trainingidx,:);
processedtest.stim = processedtest.stim(~trainingidx,:);
processedtest.spiketrain = processedtest.spiketrain(~trainingidx,:);
processedtest.stacked = processedtest.stacked(~trainingidx,:);
processedtest.intrial = processedtest.intrial(~trainingidx);

for idx = 1:nU
    sp = processed.spikes{idx};
    sptrain = sp(sp<trainmax);
    sptest = sp(sp>trainmax);
    processed.spikes{idx} = sp;
    processedtest.spikes{idx} = sptest-trainmax;
end

%Fitting stim
nS = size(processed.stim,2);
p = frames*nS;
fn_out = './monkeyresults2/run4_glm_fitting_pillow_5Hz.eps';
ggs = {};
for icell = 1:nU
    disp(num2str(icell));
    stim = processed.stim;
    stacked = processed.stacked;
    resp = processed.spikes{icell};
    sptrain = processed.spiketrain(:,icell);
    stim = stim/p;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain) - mean(stacked,1)'; 
    sta = reshape(sta,frames,[]);
    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey(sta,dt);
    gg0.tsp = resp';
    gg0.tspi = 1;
    opts = {'display', 'iter', 'maxiter', 100};
    [gg, negloglival] = MLfit_GLM_monkey(gg0,stim,opts);
    save(['./monkeyresults2/Pillow_cell_' num2str(icell) '_glm_ext_5Hz.mat'],...
        'gg'); 
    ggs{icell} = gg;
end

%Testing stim with simulation
p = frames*nS;
for icell = 1:nU
    idx = t(icell,2);
    disp(num2str(icell));
    stim = processedtest.stim;
    stim = stim/p;
    load(['./monkeyresults2/Pillow_cell_' num2str(idx) '_glm_ext_5Hz.mat']);
    Tt = size(processedtest.stim,1);
    Rt_glm = zeros(1,Tt);
    for ir = 1:nRep
        [iR_glm, vmem,Ispk] = simGLM(gg, stim);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    end
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save(['GLM_model_' num2str(idx) '_ext_5Hz.mat'],'gg','Rt_glm');
end

%Plot them
plot_filters_monkeypillow(ggs, processed, fn_out);

%Plot against my own filters 
fn_out = './monkeyresults2/run4_glm_fitting_pillow_5Hz_compare.eps';
processed_pillow = processed;
load('./monkeyresults2/run4_glm_fitting_sp_100Hz_forwardback_80pt.mat');
ggs = {};
for icell = 1:nU
    load([['./monkeyresults2/Pillow_cell_' num2str(icell) '_glm_ext_5Hz.mat']])
    ggs{icell} = gg;
end
plot_filters_monkey_compare(ggs, processed_pillow, models, data, processed, fn_out);

%Load my sims
%load('./monkeyresults2/run4_glm_fitting_sprc_100Hz_dtstm_0_01_sim.mat');
load('./monkeyresults2/run4_glm_fitting_sp_100Hz_forwardback_80pt.mat');
%Plot the spike trains... 
clf
sigma_fr = .25;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

tstart = 8001;
tend = 12000;
unitidx = 7;
unit = t(unitidx,2);
load(['GLM_model_' num2str(unit) '_ext_5Hz.mat'],'gg','Rt_glm');
tidx = tstart:tend;
truesp = processedtest.spikes(tidx,unitidx);
pillowsimsp = Rt_glm(tidx);
mysimsp = simtrains{unitidx}(tidx);
subplot(2,1,1)
hold on
plot(tidx*processed.binsize, 500*processed.grip(tidx), 'k');
plot(tidx*processed.binsize, processed.cursor(tidx,1), 'b');
plot(tidx*processed.binsize, processed.cursor(tidx,2), 'Color', [0 0.6 0]);
plot(tidx*processed.binsize, processed.cursor(tidx,3), 'r');
legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
subplot(2,1,2)
gftruesp = conv(truesp, gaussFilter_fr, 'same');
gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
gfmysimsp = conv(mysimsp, gaussFilter_fr, 'same');
plot(tidx*processed.binsize, gftruesp, tidx*processed.binsize, gfpillowsimsp, tidx*processed.binsize, gfmysimsp);
title('n rep: 20')
xlabel('seconds')
ylabel('predicted probability spiking')
legend('true spike train', 'Pillow''s GLM', 'Indep GLM')
saveplot(gcf, [fn_out '_sim'], 'eps', [6 6]);

clf
sigma_fr = .01;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

tstart = 800;
tend = 1200;
unitidx = 7;
unit = t(unitidx,2);
load(['GLM_model_' num2str(unit) '_ext_5Hz.mat'],'gg','Rt_glm');
tidx = tstart:tend;
truesp = processedtest.spikes(tidx,unitidx);
pillowsimsp = Rt_glm(tidx);
mysimsp = simtrains{unitidx}(tidx);
subplot(2,1,1)
hold on
plot(tidx*processed.binsize, 500*processed.grip(tidx), 'k');
plot(tidx*processed.binsize, processed.cursor(tidx,1), 'b');
plot(tidx*processed.binsize, processed.cursor(tidx,2), 'Color', [0 0.6 0]);
plot(tidx*processed.binsize, processed.cursor(tidx,3), 'r');
legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
subplot(2,1,2)
gftruesp = conv(truesp, gaussFilter_fr, 'same');
gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
gfmysimsp = conv(mysimsp, gaussFilter_fr, 'same');
spidx = truesp==1;
plot(tidx(spidx)*processed.binsize, truesp(spidx)-.95, '.', tidx*processed.binsize, gfpillowsimsp, tidx*processed.binsize, gfmysimsp);
title('n rep: 20')
xlabel('seconds')
ylabel('predicted probability spiking')
legend('true spike train', 'Pillow''s GLM', 'Indep GLM')
saveplot(gcf, [fn_out '_sim_zoom'], 'eps', [6 6]);