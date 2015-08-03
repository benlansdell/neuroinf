datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
nUtop = 8;
binsize = 1/100;
nK_sp = 20;
nK_stm = 6;
dt_sp = binsize;
%dt_stm = 5/100;
dt_stm = binsize;
a = 10;
const = 'on';
models = {};
nRep = 20;
fn_out = './monkeyresults/run_glm_fitting_sprc_100Hz_dtstm_0_01.eps';

processed = preprocess_monkey(datafile, binsize, 1);

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


for idx = 1:nU
	processed.unitidx = idx;
	data = filters_monkey_sprc_stm(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);
	models{idx} = MLE_glmfit(data, const);
end
plot_filters_rc_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});

%Run simulation on test stim
simtrains = {};
for idx = 1:nUtop
	idx
	processedtest.unitidx = idx;
	data = filters_monkey_sprc_stm(processedtest, nK_sp, nK_stm, a, dt_sp, dt_stm);
	sptrain = zeros(size(data.y,2), 1);
	for j = 1:nRep
		sptrain = sptrain + glmsim_rc(processedtest, models{idx}, data, 1);
	end
	sptrain = sptrain/nRep;
	simtrains{idx} = sptrain;
end
%Save the sims for comparison with other sims
save('./monkeyresults/run_glm_fitting_sprc_100Hz_dtstm_0_01_sim.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')

%Plot the spike trains... 
clf
tstart = 1;
tend = 200;
unitidx = 8;
tidx = tstart:tend;
truesp = processedtest.spikes(tidx,unitidx);
simsp = simtrains{unitidx}(tidx);
subplot(2,1,1)
plot(tidx*processed.binsize, truesp);

subplot(2,1,2)
plot(tidx*processed.binsize, truesp, '.', tidx*processed.binsize, simsp);
saveplot(gcf, [fn_out '_sim']);