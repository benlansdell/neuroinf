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
	m = models{idx};
	%m.b_hat(2:11) = m.b_hat(2:11)-0.033;
	%Simulate
	sptrain = zeros(size(testdata.y,2), 1);
	nRep = 10;
	for j = 1:nRep
		j
		[y, tspks, rho] = glmsim_rc(processed, m, testdata, 1);
		sptrain = sptrain + y;
	end
	sptrain = sptrain/nRep;
	simtrains{idx} = sptrain;

	b_hat = m.b_hat;
	%b_hat(2:11) = 0;
	%b_hat(12:end) = 20*b_hat(12:end);
	filteredstim = exp([ones(size(testdata.X,1),1), testdata.X]*b_hat');

	clf
	tstart = 10000;
	tend = 18000;
	unitidx = idx;
	tidx = tstart:tend;
	truesp = testdata.y(1,tidx);
	simsp = simtrains{unitidx}(tidx);
	filteredstim = filteredstim(tidx);
	subplot(2,1,1)
	plot(tidx*processed.binsize, 50*testdata.grip(tidx), tidx*processed.binsize, testdata.cursor(tidx,1),...
		tidx*processed.binsize, testdata.cursor(tidx,2), tidx*processed.binsize, testdata.cursor(tidx,3));	
	xlabel('s')
	legend('Grip', 'Curs x', 'Curs y', 'Curs z')
	subplot(2,1,2)
	gftruesp = conv(truesp, gaussFilter1, 'same');
	gfsimsp = conv(simsp, gaussFilter1, 'same');
	gffilteredstim = conv(filteredstim, gaussFilter1, 'same');
	plot(tidx*processed.binsize, gftruesp, '.', tidx*processed.binsize, gfsimsp, tidx*processed.binsize, gffilteredstim);
	xlabel('s')
	ylabel('Pred prob spiking')
	saveplot(gcf, [fn_out '_sim']);

	clf
	tstart = 10000;
	tend = 11000;
	unitidx = idx;
	tidx = tstart:tend;
	truesp = testdata.y(1,tidx);
	simsp = simtrains{unitidx}(tidx);
	subplot(2,1,1)
	plot(tidx*processed.binsize, 50*testdata.grip(tidx), tidx*processed.binsize, testdata.cursor(tidx,1),...
		tidx*processed.binsize, testdata.cursor(tidx,2), tidx*processed.binsize, testdata.cursor(tidx,3));	
	xlabel('s')
	legend('Grip', 'Curs x', 'Curs y', 'Curs z')
	subplot(2,1,2)
	gftruesp = conv(truesp, gaussFilter2, 'same');
	spidx = truesp > 0;
	gfsimsp = conv(simsp, gaussFilter2, 'same');
	plot(tidx(spidx)*processed.binsize, 0.1*ones(size(tidx(spidx))), 'rx', tidx*processed.binsize, gfsimsp);
	ylabel('Pred prob spiking')
	saveplot(gcf, [fn_out '_simzoom']);
end
fn_out = './monkeyresults2/run1_glm_fitting_sp_5Hz_fwdrev_trim_80pt.eps';
plot_filters_rc_monkey(models, data, processed, fn_out);

goodunits = 46-[5, 6, 7, 9, 10, 11, 12, 13, 15, 23, 25, 26, 27, 28, 29, 31, 32, 37, 38, 39, 40, 41, 42, 43];

fn_out = './monkeyresults2/run1_glm_fitting_sp_5Hz_fwdrev_trim_80pt_goodunits.eps';
plot_filters_rc_monkey_good(models, data, processed, goodunits, fn_out);

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'cursor', 'grip', 'spikes'});
%save('./monkeyresults2/run1_glm_fitting_sp_5Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed', 'simtrains', 'nRep')
save('./monkeyresults2/run1_glm_fitting_sp_5Hz_forwardback_trim_80pt.mat', 'models', 'data', 'processed')