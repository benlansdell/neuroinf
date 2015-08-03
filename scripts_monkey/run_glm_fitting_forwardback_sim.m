datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
nUtop = 8;
binsize = 1/100;
samplerate = 1/binsize;
nK_sp = 20;
nK_stm = 11;
dt_sp = binsize;
dt_stm = 200/1000;
const = 'on';
models = {};
nRep = 20;
fn_out = './monkeyresults/run_glm_fitting_sp_100Hz_fwdrev_80pt.eps';
processed = preprocess_monkey(datafile, binsize, 1);

sigma_fr = .25;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

%Only fit the top 8 units by firing rate
spikecounts = sum(processed.spikes);
tosort = zeros(nU,2);
for idx = 1:nU
	tosort(idx,:) = [spikecounts(idx),idx];
end
t=flipud(sortrows(tosort));
topunits = t(1:nUtop,2);
processed.spikes = processed.spikes(:,topunits);
processed.unitnames = processed.unitnames(topunits);

%Truncate length so fit faster
processed.cursor = processed.cursor(1:200000,:);
processed.grip = processed.grip(1:200000,:);
processed.spikes = processed.spikes(1:200000,:);

%Split into 20% test and 80% training set
processedtest = processed;
trainingidx = (1:200000)<=160000;
processed.cursor = processed.cursor(trainingidx,:);
processed.grip = processed.grip(trainingidx,:);
processed.spikes = processed.spikes(trainingidx,:);

processedtest.cursor = processedtest.cursor(~trainingidx,:);
processedtest.grip = processedtest.grip(~trainingidx,:);
processedtest.spikes = processedtest.spikes(~trainingidx,:);

%Fit
for idx = 1:nUtop
	processed.unitidx = idx;
	data = filters_monkey_sp_stmrev(processed, nK_sp, nK_stm, dt_sp, dt_stm);
	models{idx} = MLE_glmfit(data, const);
end
plot_filters_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});

%Run simulation on test stim
simtrains = {};
for idx = 1:nUtop
	idx
	processedtest.unitidx = idx;
	sptrain = zeros(size(data.y,2), 1);
	data = filters_monkey_sp_stmrev(processedtest, nK_sp, nK_stm, dt_sp, dt_stm);
	for j = 1:nRep
		sptrain = sptrain + glmsim(processedtest, models{idx}, data, 1);
	end
	sptrain = sptrain/nRep;
	simtrains{idx} = sptrain;
end
%Save the sims for comparison with other sims
save('./monkeyresults/run_glm_fitting_sp_100Hz_forwardback_80pt.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')

%Plot the spike trains... 
clf
tstart = 1;
tend = 8000;
unitidx = 4;
tidx = tstart:tend;
truesp = processedtest.spikes(tidx,unitidx);
simsp = simtrains{unitidx}(tidx);
%subplot(2,1,1)
%plot(tidx*processed.binsize, truesp);

%subplot(2,1,2)
gftruesp = conv(truesp, gaussFilter_fr, 'same');
gfsimsp = conv(simsp, gaussFilter_fr, 'same');
plot(tidx*processed.binsize, gftruesp, '.', tidx*processed.binsize, gfsimsp);
saveplot(gcf, [fn_out '_sim']);