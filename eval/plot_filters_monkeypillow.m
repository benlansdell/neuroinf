function plot_filters_monkeypillow(models, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	models = data structure output by function in ./fitting (containing fitted coefficients)
	%	fn_out = base filename to write plots to for each unit

	global RefreshRate;
	nU = 1; %number of units
	nM = length(models);
	names = {'spike history', 'curs x', 'curs y', 'curs z', 'grip'};
	nK = length(names);

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

	for i = 1:nM
		model = models{i};
		idx = 1;
		stimfilt = model.k;
		const = model.dc;
		sphist = model.ihbas*model.ih;
		for j = 1:(nK)
			%Extract data
			name = names{j};
			if j == 1
				filt = sphist(:);
				filt = flipud(filt);
			else 
				filt = stimfilt(:,j-1);
			end
			dt_filt = model.dt*RefreshRate;
			ax=axes('position',sub_pos{j,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			tt = (0:length(filt)-1)*dt_filt;
			if j == 1
				tt = (tt-max(tt))/RefreshRate/model.dt;
			end
			plot(tt, filt);
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
		str1(1) = {['Unit: ' name]};
		%str1(2) = {['Deviance: ' num2str(dev)]};
		%str1(3) = {['Degrees of freedom: ' num2str(stats.dfe)]};
		str1(2) = {['Binsize: ' num2str(processed.binsize)]};
		str1(3) = {['Seconds of training: ' num2str(size(processed.cursor,1)*processed.binsize)]};
		str1(4) = {['Number of spikes: ' num2str(sum(processed.spiketrain(:,i)))]};
		text(0.1,0.8,str1, 'FontSize', 5, 'Interpreter', 'none')
		axis off
	end
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);	
end