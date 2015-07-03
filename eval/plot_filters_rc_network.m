function plot_filters_rc_network(models, data, processed, fn_out)
	%Plot filters of a fitted network GLM model, along with other fit statistics
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
	%	data = filters_sprc_stm_network(processed, nK_sp, nK_stm);
	%	model = MLE_glmfit(data, const);
	%	plot_filters_rc_network(models, data, processed, fn_out);

	nM = length(models);
	nK_stm = data.nK_stm;
	nU = size(data.k,1)-1;

	%Test run
	%nM = 3;

	clf;
	maxbstm = 0;
	minbstm = 0;
	maxbsp = 0;
	minbsp = 0;
	%Same color scale for all
	for idx = 1:nM
		model = models{idx};
		k = data.k; %filter names and indices
		b_hat = model.b_hat(1,:);
		for j = 1:nU
			spidx = data.k{j,2};
			b_sp = data.spbasis*b_hat(1,1+spidx)';
			if model.converged & model.conditioned 
				maxbsp = max(maxbsp, max(b_sp));
				minbsp = min(minbsp, min(b_sp));
			end
		end
		stmidx = k{nU+1,2};
		b_stm = b_hat(1,1+stmidx);
		if model.converged & model.conditioned 
			maxbstm = max(maxbstm, max(b_stm));
			minbstm = min(minbstm, min(b_stm));
		end
	end

	%Prepare subplots
	plotheight=nM+3;
	plotwidth=50;
	subplotsx=nU+nK_stm;
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
		stmidx = k{nU+1,2};
		b_stm = b_hat(1,1+stmidx);
		nX = length(stmidx)/nK_stm;
		Nv = processed.Nv;
		nP = nU + nK_stm;
		for j = 1:nU
			ax=axes('position',sub_pos{j,idx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			spidx = data.k{j,2};
			b_sphist = data.spbasis*b_hat(1,1+spidx)';
			if model.conditioned & model.converged
				plot(b_sphist)
			end
			ylim([minbsp, maxbsp])
			axis off
		end
		%ylabel('spike history')
		colormap(pink);
		for j = 1:nK_stm
			ax=axes('position',sub_pos{nU+j,idx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			frame = reshape(b_stm(nX*(j-1)+1:nX*j), Nv, []);
			if model.conditioned & model.converged
				h = pcolor(frame);
			end
			set(h, 'EdgeColor', 'none');
			axis off
			caxis([minbstm, maxbstm])
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
		ax=axes('position',sub_pos{nU+idx,nM+1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		axis off
		if idx == 1
			str1(1) = {'Stimulus'};
			str1(2) = {['t-' num2str(nK_stm+1-idx) '\Delta t']};
		else
			str1(1) = {['t-' num2str(nK_stm+1-idx) '\Delta t']};
		end
		text(0.1,0.8,str1)
	end
	caxis([minbstm, maxbstm])
	colorbar

	%save eps
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);

end