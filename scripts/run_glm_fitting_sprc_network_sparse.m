responsefiles = struct2cell(dir('./data/*.isk'));
responsefiles = responsefiles(1,:);
stimfile = './data/whitenoise.raw';
fn_out = './run_glm_fitting_sprc_network_sparse.eps';

nK_sp = 20;
nK_stm = 6;
const = 'on';
binsize = 1;
a = 11;
%As determined by optimallambda.m
lambda = 30;

models = {};
for idx = 1:length(responsefiles)
	display(['Fitting unit ' num2str(idx)])
	processed = preprocess(stimfile, responsefiles, binsize, idx);
	data = filters_sprc_stm_network(processed, nK_sp, nK_stm, a);
	models{idx} = MLE_SD_rc_reg(data, lambda, const);
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_glm_fitting_sprc_network_sparse.mat', 'models', 'data', 'processed')
plot_filters_rc_network(models, data, processed, fn_out);