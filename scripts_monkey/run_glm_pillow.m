global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100 ; 
dt = 0.1;
datafile = './data/mabel_reaching_5-4-10.mat';
%Test run...
%nU = 4;
nU = 45;
nS = 4;
frames = 6;
binsize = 1/RefreshRate;
processed = preprocess_monkey_pillow(datafile, binsize, dt, frames);    
nS = size(processed.stim,2);
p = frames*nS;
fn_out = './monkeyresults/run_glm_fitting_pillow.eps';
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
    save(['./monkeyresults/Pillow_cell_' num2str(icell) '_glm.mat'],...
        'gg'); 
    ggs{icell} = gg;
end

%Plot them
plot_filters_monkeypillow(ggs, processed, fn_out);

%Plot against my own filters 
fn_out = './monkeyresults/run_glm_fitting_pillow_compare_dtstm_0_01.eps';
processed_pillow = processed;
load('./monkeyresults/run_glm_fitting_sprc_100Hz_dtstm_0_01.mat');
ggs = {};
for icell = 1:nU
    load([['./monkeyresults/Pillow_cell_' num2str(icell) '_glm.mat']])
    ggs{icell} = gg;
end
plot_filters_monkey_compare(ggs, processed_pillow, models, data, processed, fn_out);
