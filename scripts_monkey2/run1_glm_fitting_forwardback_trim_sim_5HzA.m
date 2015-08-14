datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
samplerate = 1/binsize;
nK_sp = 40;
nK_stm = 5;
a = 10;
dt_sp = binsize;
dt_stm = 200/1000;
const = 'on';
percentage = 80;
models = {};
nRep = 20;
fn_out = './monkeyresults2/run1_glm_fitting_sp_5Hz_fwdrev_trim_80pt.eps';
processed = preprocess_monkey(datafile, binsize, 1);
nB = size(processed.cursor,1);


sigma_fr = .25;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter1 = exp(-x.^2/(2*sigma_fr^2));
gaussFilter1 = gaussFilter1/sum(gaussFilter1);

sigma_fr = .02;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter2 = exp(-x.^2/(2*sigma_fr^2));
gaussFilter2 = gaussFilter2/sum(gaussFilter2);

%Fit and simulate
for idx = 1:nU
	fn_out = ['./monkeyresults2/run1_glm_fitting_sp_5Hz_fwdrev_trim_80pt_unit_' num2str(idx) '.eps'];
	%Test run
	%idx = 4;
	processed.unitidx = idx;
	data = filters_monkey_sprc_stmrev_trim(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);
	[data, testdata] = splitforCV(data, percentage);
	models{idx} = MLE_glmfit(data, const);
end
fn_out = './monkeyresults2/run1_glm_fitting_sp_5Hz_fwdrev_trim_80pt.eps';
plot_filters_rc_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
%save('./monkeyresults2/run1_glm_fitting_sp_5Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')
save('./monkeyresults2/run1_glm_fitting_sp_5Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed')