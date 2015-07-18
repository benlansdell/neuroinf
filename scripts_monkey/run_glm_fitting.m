datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
nK_sp = 20;
nK_stm = 6;
dt_sp = binsize;
dt_stm = 5/100;
const = 'on';
models = {};
fn_out = './monkeyresults/run_glm_fitting_sp.eps';

processed = preprocess_monkey(datafile, binsize, 1);
for idx = 1:nU
	processed.unitidx = idx;
	data = filters_monkey_sp_stm(processed, nK_sp, nK_stm, dt_sp, dt_stm);
	models{idx} = MLE_glmfit(data, const);
end
plot_filters_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
save('./monkeyresults/run_glm_fitting_sp.mat', 'models', 'data', 'processed')
