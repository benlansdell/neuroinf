datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
samplerate = 1/binsize;
nK_sp = 40;
nK_stm = 10;
a = 10;
dt_sp = binsize;
dt_stm = 100/1000;
const = 'on';
percentage = 80;
models = {};
nRep = 20;
fn_out = './monkeyresults2/run5_glm_fitting_sp_10Hz_fwdrev_trim_80pt.eps';
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
	%Test run
	%idx = 7;
	processed.unitidx = idx;
	data = filters_monkey_sprc_stmrev_trim(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);
	[data, testdata] = splitforCV(data, percentage);
	models{idx} = MLE_glmfit(data, const);
	%models{idx}.b_hat(2:11) = models{idx}.b_hat(2:11)-0.05;
	%Simulate
	%sptrain = zeros(size(testdata.y,2), 1);
	%nRep = 10;
	%for j = 1:nRep
	%	j
	%	sptrain = sptrain + glmsim_rc(processed, models{idx}, testdata, 1);
	%end
	%sptrain = sptrain/nRep;
	%simtrains{idx} = sptrain;

	%clf
	%tstart = 1;
	%tend = 8000;
	%unitidx = idx;
	%tidx = tstart:tend;
	%truesp = testdata.y(1,tidx);
	%simsp = simtrains{unitidx}(tidx);
	%subplot(2,1,1)
	%plot(tidx*processed.binsize, 50*testdata.grip(tidx), tidx*processed.binsize, testdata.cursor(tidx,1),...
	%	tidx*processed.binsize, testdata.cursor(tidx,2), tidx*processed.binsize, testdata.cursor(tidx,3));	
	%xlabel('s')
	%legend('Grip', 'Curs x', 'Curs y', 'Curs z')
	%subplot(2,1,2)
	%gftruesp = conv(truesp, gaussFilter1, 'same');
	%gfsimsp = conv(simsp, gaussFilter1, 'same');
	%plot(tidx*processed.binsize, gftruesp, '.', tidx*processed.binsize, gfsimsp);
	%xlabel('s')
	%ylabel('Pred prob spiking')
	%saveplot(gcf, [fn_out '_sim']);
%
	%clf
	%tstart = 2000;
	%tend = 3000;
	%unitidx = idx;
	%tidx = tstart:tend;
	%truesp = testdata.y(1,tidx);
	%simsp = simtrains{unitidx}(tidx);
	%subplot(2,1,1)
	%plot(tidx*processed.binsize, 50*testdata.grip(tidx), tidx*processed.binsize, testdata.cursor(tidx,1),...
	%	tidx*processed.binsize, testdata.cursor(tidx,2), tidx*processed.binsize, testdata.cursor(tidx,3));	
	%xlabel('s')
	%legend('Grip', 'Curs x', 'Curs y', 'Curs z')
	%subplot(2,1,2)
	%gftruesp = conv(truesp, gaussFilter2, 'same');
	%spidx = truesp > 0;
	%gfsimsp = conv(simsp, gaussFilter2, 'same');
	%plot(tidx(spidx)*processed.binsize, 0.1, 'rx', tidx*processed.binsize, gfsimsp);
	%ylabel('Pred prob spiking')
	%saveplot(gcf, [fn_out '_simzoom']);
end
plot_filters_rc_monkey(models, data, processed, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
%save('./monkeyresults2/run5_glm_fitting_sp_10Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')
save('./monkeyresults2/run5_glm_fitting_sp_10Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed')