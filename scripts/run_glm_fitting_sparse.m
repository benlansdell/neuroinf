responsefiles = struct2cell(dir('./data/*.isk'));
responsefiles = responsefiles(1,:);
responsefiles = sortnumerical(responsefiles);
stimfile = './data/whitenoise.raw';
fn_out = './run_glm_fitting_sp_sparse.eps';

nK_sp = 6;
nK_stm = 6;
const = 'on';
binsize = 1;
models = {};
lambda = 35;

for idx = 1:length(responsefiles)
	display(['Fitting unit ' num2str(idx)])
	rf = responsefiles(idx);
	processed = preprocess(stimfile, rf, binsize, 1);
	data = filters_sp_stm(processed, nK_sp, nK_stm);
	models{idx} = MLE_SD_reg(data, lambda, const);
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_glm_fitting_sp_sparse.mat', 'models', 'data', 'processed')
plot_filters(models, data, processed, fn_out);