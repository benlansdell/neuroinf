datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
nK_sp = 20;
nK_stm = 6;
dt_sp = binsize;
dt_stm = 5/100;
a = 10;
const = 'on';
models = {};
fn_out = './monkeyresults/run_glm_fitting_sprc_network.eps';


global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100 ; 
dt = 0.1;
nS = 4;
frames = 6;
processedp = preprocess_monkey_pillow(datafile, 1/RefreshRate, dt, frames);   

processed = preprocess_monkey(datafile, binsize, 1);
for idx = 1:nU
	processed.unitidx = idx;
	data = filters_monkey_sprc_stm_network(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);
	models{idx} = MLE_glmfit(data, const);
end

load('./monkeyresults/run_glm_fitting_sprc_network.mat')
plot_filters_rc_network_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
save('./monkeyresults/run_glm_fitting_sprc_network.mat', 'models', 'data', 'processed')


%Compare to Pillow's network
modelsp = {};
for icell = 1:nU
    load(['./monkeyresults/Pillow_cell_' num2str(icell) '_glm.mat']); 
    modelsp{icell} = gg;
end

plot_filters_rc_network_monkey_compare(models, data, processed, modelsp, processedp, [fn_out '_compare']);
