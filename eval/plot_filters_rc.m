function plot_filters_rc(models, data, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%	fn_out = base filename to write plots to for each unit
	%
	%Test code:
	%
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 20;
	%	nK_stm = 6;
	%	const = 'on';
	%	fn_out = './testfilters.eps';
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm(processed, nK_sp, nK_stm);
	%	model = MLE_glmfit(data, const);
	%	plot_filters_rc(model, data, processed, fn_out);

	nM = length(models);
	%Test run
	%nM = 3;
	clf;
	maxbstm = 0;
	minbstm = 0;
	%Same color scale for all
	for idx = 1:nM
		model = models{idx};
		k = data.k; %filter names and indices
		b_hat = model.b_hat(1,:);
		stmidx = k{2,2};
		b_stm = b_hat(1,1+stmidx);
		maxbstm = max(maxbstm, max(b_stm));
		minbstm = min(minbstm, min(b_stm));
	end

	for idx = 1:nM
		idx
		model = models{idx};
		k = data.k; %filter names and indices
		b_hat = model.b_hat(1,:);
		dev = model.dev{1};
		stats = model.stats{1};
		const = b_hat(1);
		spidx = data.sp_hist;
		stmidx = k{2,2};
		b_sphist = data.spbasis*b_hat(1,1+spidx)';
		b_stm = b_hat(1,1+stmidx);
		nK_stm = data.nK_stm;
		nX = length(stmidx)/nK_stm;
		Nv = processed.Nv;
		nP = 1 + nK_stm;
		subplot(nM+1,nP,(nP*(idx-1)+1))
		if model.converged & model.conditioned
			plot(b_sphist)
		else
			plot(b_sphist, '--')
		end
		%ylabel('spike history')
		colormap(pink);
		for j = 1:nK_stm
			subplot(nM+1,nP,(nP*(idx-1)+1+j))
			frame = reshape(b_stm(nX*(j-1)+1:nX*j), Nv, []);
			h = pcolor(frame);
			set(h, 'EdgeColor', 'none');
			axis off
			caxis([min(b_stm), max(b_stm)])
		end
	end
	subplot(nM+1, nP, nP*nM+1)
	axis off
	str1(1) = {['Spike history']};
	str1(2) = {['Solid=conv,cond']};
	str1(3) = {['Dashed=not']};
	text(0,0.8,str1)

	for idx = 1:nK_stm
		str1 = {};
		subplot(nM+1, nP, nP*nM+1+idx)
		axis off
		if idx == 1
			str1(1) = {'Stimulus'};
			str1(2) = {['t-' num2str(nK_stm+1-idx) '\Delta t']};
		else
			str1(1) = {['t-' num2str(nK_stm+1-idx) '\Delta t']};
		end
		text(0.1,0.8,str1)
	end
	%caxis([minbstm, maxbstm])
	%colorbar

	%save eps
	saveplot(gcf, fn_out, 'eps', [10,30]);
end



