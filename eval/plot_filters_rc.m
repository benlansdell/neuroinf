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
	nK_stm = data.nK_stm;
	nU = 1;

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

	%Prepare subplots
	plotheight=nM+3;
	plotwidth=8;
	subplotsx=nU+nK_stm+1;
	subplotsy=nM+1;   
	leftedge=1.2;
	rightedge=0.4;   
	topedge=1;
	bottomedge=1.5;
	spacex=0.2;
	spacey=0.2;
	fontsize=5;    
	sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

	%setting the Matlab figure
	f=figure;
	clf(f);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [plotwidth plotheight]);
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

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
		%subplot(nM+1,nP,(nP*(idx-1)+1))
		ax=axes('position',sub_pos{1,nM+1-idx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		axis off
		str1(1) = {['Unit ' num2str(idx)]};
		str1(2) = {['# spikes:']};
		str1(3) = {[num2str(model.nspikes)]};
		text(0,0.8,str1)

		ax=axes('position',sub_pos{2,nM+1-idx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

		if model.converged & model.conditioned
			plot(processed.binsize*((1:length(b_sphist))-length(b_sphist)-1), b_sphist)
		else
			plot(processed.binsize*((1:length(b_sphist))-length(b_sphist)-1), b_sphist, '--')
		end
		%ylabel('spike history')
		colormap(pink);
		for j = 1:nK_stm
			ax=axes('position',sub_pos{2+j,nM+1-idx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			frame = reshape(b_stm(nX*(j-1)+1:nX*j), Nv, []);
			h = pcolor(frame);
			set(h, 'EdgeColor', 'none');
			axis off
			caxis([min(b_stm), max(b_stm)])
		end
	end
	ax=axes('position',sub_pos{1,nM+1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
	axis off
	str1(1) = {['Spike history']};
	str1(2) = {['Solid=conv,cond']};
	str1(3) = {['Dashed=not']};
	text(0,0.8,str1)

	for idx = 1:nK_stm
		str1 = {};
		ax=axes('position',sub_pos{2+idx,nM+1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
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
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);
end



