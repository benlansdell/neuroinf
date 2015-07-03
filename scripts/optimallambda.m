%Find optimal value of lambda according to deviance on a test dataset
%responsefiles = struct2cell(dir('./data/*.isk'));
%responsefiles = responsefiles(1,1:6);
responsefiles = {'whitenoisec1.isk', 'whitenoisec2.isk', 'whitenoisec3.isk',...
	'whitenoisec4.isk','whitenoisec5.isk','whitenoisec6.isk'};
stimfile = './data/whitenoise.raw';
fn_out = './optimallambda.eps';

nK_sp = 20;
nK_stm = 6;
const = 'on';
binsize = 1;
a = 10;
lambdas = [0 5 10 15 20 25 30 35 40 45 50];

models = {};
devs = [];
for l = 1:length(lambdas)
	lambda = lambdas(l)
	for idx = 1:length(responsefiles)
		display(['Fitting unit ' num2str(idx)])
		processed = preprocess(stimfile, responsefiles, binsize, idx);
		%Truncate to a shorter time
		training = processed;
		training.spikes = training.spikes(1:35000,:);
		training.stim = training.stim(1:35000,:);
		%Fit model to training
		data = filters_sprc_stm_network(training, nK_sp, nK_stm, a);
		models{l,idx} = MLE_SD_rc_reg(data, lambda, const);
		%Find deviance on test dataset
		processed.spikes = processed.spikes(35001:70000,:);
		processed.stim = processed.stim(35001:70000,:);
		data = filters_sprc_stm_network(processed, nK_sp, nK_stm, a);
		devs(l,idx) = deviance_network(models{l,idx}, data);
		devs(l,idx)
	end
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_optimallambda.mat', 'models', 'data', 'processed')

%Plot results