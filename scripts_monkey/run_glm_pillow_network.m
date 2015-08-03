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
fn_out = './monkeyresults/run_glm_fitting_pillow_network.eps';
ggs = {};
for icell = 1:nU
    disp(num2str(icell));
    stim = processed.stim;
    stacked = processed.stacked;
    resp = processed.spikes{icell};
    nicell = [(1:icell-1), (icell+1:nU)];
    %Add coupling to the other spike trains
    coupled = processed.spikes(nicell);
    for idx = 1:length(coupled)
    	coupled{idx} = coupled{idx}';
    end
    sptrain = processed.spiketrain(:,icell);
    stim = stim/p;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain) - mean(stacked,1)'; 
    sta = reshape(sta,frames,[]);
    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey(sta,dt);
    gg0.tsp = resp';
    gg0.tspi = 1;
    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);
    opts = {'display', 'iter', 'maxiter', 100};
    [gg, negloglival] = MLfit_GLM_monkey(gg0,stim,opts);
    save(['./monkeyresults/Pillow_cell_' num2str(icell) '_glm_network.mat'],...
        'gg'); 
    ggs{icell} = gg;
end
plot_filters_monkeypillownetwork(ggs, processed, fn_out);