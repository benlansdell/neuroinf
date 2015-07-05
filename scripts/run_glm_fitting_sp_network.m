responsefiles = struct2cell(dir('./data/*.isk'));
responsefiles = responsefiles(1,:);
responsefiles = sortnumerical(responsefiles);
stimfile = './data/whitenoise.raw';
fn_out = './run_glm_fitting_sp_network.eps';

nK_sp = 6;
nK_stm = 6;
const = 'on';
binsize = 1;
models = {};

for idx = 1:length(responsefiles)
	display(['Fitting unit ' num2str(idx)])
	processed = preprocess(stimfile, responsefiles, binsize, idx);
	data = filters_sp_stm_network(processed, nK_sp, nK_stm);
	models{idx} = MLE_glmfit(data, const);
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_glm_fitting_sp_network.mat', 'models', 'data', 'processed')
plot_filters_network(models, data, processed, fn_out);