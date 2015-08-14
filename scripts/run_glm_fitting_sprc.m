responsefiles = struct2cell(dir('./data/*.isk'));
responsefiles = responsefiles(1,:);
responsefiles = sortnumerical(responsefiles);
stimfile = './data/whitenoise.raw';
fn_out = './run_glm_fitting_sprc_a_10.eps';

nK_sp = 6;
nK_stm = 6;
const = 'on';
binsize = 1/30;
models = {};
a = 10;

for idx = 1:length(responsefiles)
	display(['Fitting unit ' num2str(idx)])
	rf = responsefiles(idx);
	processed = preprocess(stimfile, rf, binsize, 1);
	data = filters_sprc_stm(processed, nK_sp, nK_stm, a);
	models{idx} = MLE_glmfit(data, const);
end

%Save models, and other data
data = rmfield(data, {'X', 'y'});
processed = rmfield(processed, {'stim', 'spikes'});
save('./run_glm_fitting_sprc.mat', 'models', 'data', 'processed')
plot_filters_rc(models, data, processed, fn_out);

%Simulate models
nRep = 20;
for idx = 1:length(responsefiles)
	display(['Simulating unit ' num2str(idx)])
	rf = responsefiles(idx);
	processed = preprocess(stimfile, rf, binsize, 1);
	data = filters_sprc_stm(processed, nK_sp, nK_stm, a);
	sptrain = zeros(size(data.y))';
	load('./run_glm_fitting_sprc.mat', 'models')
	for j = 1:nRep
		j
		model = models{j};
		sptrain = sptrain + glmsim(processed, model, data, 1);
	end
	sptrains{idx} = sptrain/nRep;
end
save('./run_glm_fitting_sprc_sims.mat', 'models', 'sptrains');

sigma = .01;
sigma = sigma/binsize;
sz = sigma*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter = exp(-x.^2/(2*sigma^2));
gaussFilter = gaussFilter/sum(gaussFilter);

%Simulate Pillow's code for comparison:
pX = 100;
pT = 6;
mses_matlab = [];
mses_pillow = [];
for idx = 1:length(responsefiles)
	load(['./data/Retina_cell_' num2str(idx) '_glmamended_long.mat'])
	load(['./data/Retina_cell_' num2str(idx) '_stim_resp_long.mat'])
	p = pX*pT;
	St = St/p;
	S1t = St(:,1:pX);
	%Rt_glm = zeros(1,Tt);
	%for ir = 1:nRep
	%    [iR_glm, vmem,Ispk] = simGLM(gg, S1t);
	%    Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
	%end
	%Rt_glm = Rt_glm'/nRep + 1e-8;
	ii = 1:200;
	ii_m = ii+floor(size(data.y,2)*.8);
	tt = ii*processed.binsize;
	plot(tt, sptrains{idx}(ii_m), tt, 0.1*(Rt(ii)>0), '.', tt, Rt_glm(ii));
	legend('MATLAB', 'true', 'Pillow')
	saveplot(gcf, ['./simplots/compare_pillow_matlab_glm_cell_' num2str(idx) '.eps'])

	%mses_matlab(idx) = mean((sptrains{idx}(ii_m(1):end)-Rt).^2);
	%mses_pillow(idx) = mean((Rt_glm-Rt).^2);
end