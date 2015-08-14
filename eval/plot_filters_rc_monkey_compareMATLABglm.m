function plot_filters_rc_monkey_compareMATLABglm(mdls, data, processed, mdls2, data2, processed2, fn_out)

	%nU = size(processed.cursor,1); %number of units
	nU = length(mdls); %number of units
	k = data.k; %filter names and indices
	nK = size(k,1); %number of filters
	nP = 3; %number of things to plot about each fitted filter
	%For each unit, plot filters fit
	h = figure;
	if ~iscell(mdls)
		models{1} = mdls;
		models2{1} = mdls2;
	else
		models = mdls;
		models2 = mdls2;
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

	%Before plotting, collect the largest cursor and grip values to set scale by
	ymincurs = 0;
	ymaxcurs = 0;
	ymingrip = 0;
	ymaxgrip = 0;
	for i = 1:nM
		model = models{i};
		model2 = models2{i};
		idx = 1;
		b_hat = model.b_hat(idx,:);
		b_hat2 = model2.b_hat(idx,:);
		if model.conditioned == 0
			continue
		end
		if model.converged == 0
			continue
		end
		for j = 1:nK
			%Extract data
			name = k{j,1};
			filt = b_hat(k{j,2}+1);
			filt2 = b_hat2(k{j,2}+1);
			if j > 1 & j < 5
				ymincurs = min(ymincurs, min(filt));
				ymaxcurs = max(ymaxcurs, max(filt));
				ymincurs = min(ymincurs, min(filt2));
				ymaxcurs = max(ymaxcurs, max(filt2));
			elseif j ==5
				ymingrip = min(ymingrip, min(filt));
				ymaxgrip = max(ymaxgrip, max(filt));
				ymingrip = min(ymingrip, min(filt2));
				ymaxgrip = max(ymaxgrip, max(filt2));
			end
		end
	end

	clf;
	for i = 1:nM
		model = models{i};
		idx = 1;
		b_hat = model.b_hat(idx,:);
		dev = model.dev{idx};
		stats = model.stats{idx};
		const = b_hat(1);
		if model.conditioned == 0
			continue
		end
		if model.converged == 0
			continue
		end		

		model2 = models2{i};
		b_hat2 = model2.b_hat(idx,:);
		dev2 = model2.dev{idx};
		stats2 = model2.stats{idx};
		const2 = b_hat2(1);
		if model2.conditioned == 0
			continue
		end
		if model2.converged == 0
			continue
		end		

		for j = 1:nK
			%Extract data
			name = k{j,1};
			filt = b_hat(k{j,2}+1);
			filt2 = b_hat2(k{j,2}+1);
			%Raised cosine basis functions... back to temporal basis
			if j ==1 
				filt = data.spbasis*filt';
				filt2 = data.spbasis*filt2';
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
				ymin = min([min(filt), min(filt2)])*1.2;
				ymax = max([max(filt), max(filt2)])*1.2;
			elseif j > 1 & j < 5
				tt = tt-max(tt)/2;
				ymin = ymincurs*1.2;
				ymax = ymaxcurs*1.2;
			elseif j == 5
				tt = tt-max(tt)/2;
				ymin = ymingrip*1.2;
				ymax = ymaxgrip*1.2;
			end
			hold on
			%ymin = min(filt-se)*1.2;
			%ymax = max(filt+se)*1.2;
			plot(tt, filt, tt, filt2, '--');
			ylim([ymin ymax]);
			if length(tt) > 1
				xlim([min(tt) max(tt)]);
			end
			if i == nM
				title(name);
			end
			if i == 1
				xlabel('time (ms)');
			end
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
			if isfield(data, 'y')
				str1(6) = {['Seconds of training: ' num2str(size(data.y,2)*processed.binsize)]};
			end
			if isfield(processed, 'spikes')
				str1(7) = {['Number of spikes: ' num2str(sum(processed.spikes(:,i)))]};
			end
			text(0.1,0.8,str1, 'FontSize', 5, 'Interpreter', 'none')
			axis off
		end
	end
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);	
end