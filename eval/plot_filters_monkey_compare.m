function plot_filters_monkey_compare(pillow_models, pillow_processed, mdls, data, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics

	global RefreshRate
	nU = size(processed.unitnames,1); %number of units
	k = data.k; %filter names and indices
	nK = size(k,1); %number of filters
	nP = 3; %number of things to plot about each fitted filter
	%For each unit, plot filters fit
	h = figure;
	if ~iscell(mdls)
		models{1} = mdls;
	else
		models = mdls;
	end
	nM = length(models);


	%Prepare subplots
	plotheight=nM+4;
	plotwidth=10;
	subplotsx=6;
	subplotsy=nM+1;   
	leftedge=1.2;
	rightedge=0.4;   
	topedge=1;
	bottomedge=1.5;
	spacex=0.2;
	spacey=0.4;
	fontsize=5;    
	sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

	%setting the Matlab figure
	f=figure;
	clf(f);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [plotwidth plotheight]);
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

	clf;
	for i = 1:nM
		model = models{i};

		idx = 1;
		%Indep filters
		b_hat = model.b_hat(idx,:);
		dev = model.dev{idx};
		stats = model.stats{idx};
		const = b_hat(1);

		%Pillow's filters
		pmodel = pillow_models{i};
		pstimfilt = pmodel.k;
		pconst = pmodel.dc;
		psphist = pmodel.ihbas*pmodel.ih;

		for j = 1:nK
			%Extract data
			name = k{j,1};
			filt = b_hat(k{j,2}+1);
			%If available, find statistics of fit, otherwise set these to zero
			if isfield(stats, 'se')	
				se = stats.se(k{j,2}+1)';
				tstat = stats.t(k{j,2}+1);
				pval = stats.p(k{j,2}+1);
			else
				se = zeros(size(k{j,2}));
				tstat = zeros(size(k{j,2}));
				pval = zeros(size(k{j,2}));
			end
			if j ==1 
				filt = data.spbasis*filt';
				se = zeros(size(filt));
			end
			dt_filt = k{j,3};
			%If filter length is zero skip this one
			if length(k{j,2}) < 1
				continue
			end

			%Plot filter plus/minus SE
			%subplot(nP, nK+1, j)
			ax=axes('position',sub_pos{j,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			tt = (0:length(filt)-1)*dt_filt*1000;
			if j == 1
				tt = tt-max(tt);
			end
			hold on
			ymin = min(filt-se)*1.2;
			ymax = max(filt+se)*1.2;
			area(tt, filt+se, ymin, 'FaceColor', [0.8 0.8 0.8])
			area(tt, filt-se, ymin, 'FaceColor', [1 1 1])
			plot(tt, filt);
			%ylim([ymin ymax]);
			if length(tt) > 1
				xlim([min(tt) max(tt)]);
			end
			if i == nM
				title(name);
			end
			if i == 1
				xlabel('time (ms)');
			end

			%Add pillows filters to plot:
			%Spike history filter
			if j == 1
				pfilt = psphist(:);
				pfilt = flipud(pfilt);
			else
				%Stim filters
				pfilt = pstimfilt(:,j-1);
			end
			tt = (0:length(pfilt)-1);
			if j > 1
				tt = tt/pmodel.dt;
				%Normalize
				pfilt = pfilt*pmodel.dt;
			else
				tt = tt-max(tt);
			end
			plot(tt, pfilt, 'r');
		end

		if ischar(processed.unitnames)
			name = processed.unitnames;
		else
			name = processed.unitnames{i};
		end

		%Plot information about each subplot
		%subplot(nP, nK+1, (nK+1))
		ax=axes('position',sub_pos{nK+1,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		if isfield(stats, 'dfe')
			str1(1) = {['Unit: ' name]};
			str1(2) = {['Deviance: ' num2str(dev)]};
			str1(3) = {['Degrees of freedom: ' num2str(stats.dfe)]};
			str1(4) = {['Estimated dispersion: ' num2str(stats.sfit)]};
			str1(5) = {['Binsize: ' num2str(processed.binsize)]};
			str1(6) = {['Seconds of training: ' num2str(size(data.cursor,1)*processed.binsize)]};
			str1(7) = {['Number of spikes: ' num2str(sum(pillow_processed.spiketrain(:,i)))]};
			text(0.1,0.8,str1, 'FontSize', 5, 'Interpreter', 'none')
			axis off
		end
	end
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);	
end