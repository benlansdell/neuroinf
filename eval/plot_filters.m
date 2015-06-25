function plot_filters(model, data, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	model = data structure output by function in ./fitting (containing fitted coefficients)
	%	data = data structure output by ./models containing data used for fit
	%	processed = data structure output by ./preprocess containing processed raw data
	%	fn_out = base filename to write plots to for each unit
	%

	k = data.k; %filter names and indices
	b_hat = model.b_hat(1,:);
	dev = model.dev{1};
	stats = model.stats{1};
	const = b_hat(1);
	spidx = data.sp_hist;
	stmidx = k{2,2};
	b_sphist = b_hat(1,1+spidx);
	b_stm = b_hat(1,1+stmidx);
	maxbstm = max(b_stm);
	minbstm = min(b_stm);
	nK_stm = data.nK_stm;
	nX = length(stmidx)/nK_stm;
	Nv = processed.Nv;
	nP = 2 + nK_stm;
	clf;
	subplot(1,nP,1)
	plot(b_sphist)
	xlabel('bin')
	ylabel('spike history')
	colormap(pink);
	for idx = 1:nK_stm
		subplot(1,nP,1+idx)
		frame = reshape(b_stm(nX*(idx-1)+1:nX*idx), Nv, []);
		pcolor(frame)
		caxis([minbstm, maxbstm])
	end
	subplot(1,nP,nP);
	caxis([minbstm, maxbstm])
	axis off;
	%ylabel('asdf')
	colorbar('westoutside');

	%save eps
	saveplot(gcf, fn_out, 'eps', [12,6]);
end