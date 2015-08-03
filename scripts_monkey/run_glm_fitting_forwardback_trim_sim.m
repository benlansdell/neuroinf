datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
samplerate = 1/binsize;
nK_sp = 40;
nK_stm = 6;
a = 10;
dt_sp = binsize;
dt_stm = 100/1000;
const = 'on';
percentage = 80;
models = {};
nRep = 20;
fn_out = './monkeyresults/run_glm_fitting_sp_100Hz_fwdrev_trim_80pt.eps';
processed = preprocess_monkey(datafile, binsize, 1);
nB = size(processed.cursor,1);

%Fit and simulate
for idx = 1:nU
	processed.unitidx = idx;
	data = filters_monkey_sprc_stmrev_trim(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);
	[data, testdata] = splitforCV(data, percentage);
	models{idx} = MLE_glmfit(data, const);
	%Simulate
	%sptrain = zeros(size(testdata.y,2), 1);
	%for j = 1:nRep
	%	sptrain = sptrain + glmsim(processed, models{idx}, testdata, 1);
	%end
	%sptrain = sptrain/nRep;
	%simtrains{idx} = sptrain;
end
plot_filters_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
save('./monkeyresults/run_glm_fitting_sp_100Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')

%Plot the simulated spike trains... 
%sigma_fr = .25;
%sigma_fr = sigma_fr*samplerate;
%sz = sigma_fr*3*2;
%x = linspace(-sz/2, sz/2, sz);
%gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
%gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

%clf
%tstart = 1;
%tend = 8000;
%unitidx = 4;
%tidx = tstart:tend;
%truesp = processedtest.spikes(tidx,unitidx);
%simsp = simtrains{unitidx}(tidx);
%%subplot(2,1,1)
%%plot(tidx*processed.binsize, truesp);
%
%%subplot(2,1,2)
%gftruesp = conv(truesp, gaussFilter_fr, 'same');
%gfsimsp = conv(simsp, gaussFilter_fr, 'same');
%plot(tidx*processed.binsize, gftruesp, '.', tidx*processed.binsize, gfsimsp);
%saveplot(gcf, [fn_out '_sim']);