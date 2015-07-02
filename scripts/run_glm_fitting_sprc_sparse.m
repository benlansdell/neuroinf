responsefiles = struct2cell(dir('./data/*.isk'));
responsefiles = responsefiles(1,:);
stimfile = './data/whitenoise.raw';
fn_out = './run_glm_fitting_sprc_sparse.eps';

nK_sp = 20;
nK_stm = 6;
const = 'on';
binsize = 1;
models = {};
a = 10;
lambda = 30;

for idx = 1:length(responsefiles)
	display(['Fitting unit ' num2str(idx)])
	%Find the idx of the idxth file
	unit = find(cellfun(@(x) ~isempty(findstr(x, ['sec' num2str(idx) '.isk'])),...
	 responsefiles));
	rf = responsefiles(unit);
	processed = preprocess(stimfile, rf, binsize, idx);
	data = filters_sprc_stm(processed, nK_sp, nK_stm, a);
	models{idx} = MLE_SD_rc_reg(data, lambda, const);
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_glm_fitting_sprc_sparse.mat', 'models', 'data', 'processed')
plot_filters_rc(models, data, processed, fn_out);