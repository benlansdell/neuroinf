function plot_filters_monkeypillownetwork(models, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	models = data structure output by function in ./fitting (containing fitted coefficients)
	%	fn_out = base filename to write plots to for each unit

	global RefreshRate;
	nU = 1; %number of units
	nM = length(models);
	names = {'spike history', 'curs x', 'curs y', 'curs z', 'grip'};

	%Only plot the top 10 units by firing rate:
	nUtop = 8;
	nK = nUtop+5;
	topunits = [];
	for idx = 1:nM
		topunits(idx,2) = idx;
		topunits(idx,1) = sum(processed.spiketrain(:,idx));
	end
	topunits = flipud(sortrows(topunits));
	topunits = topunits(1:nUtop,2);

	%Prepare subplots
	plotheight=nUtop+5;
	plotwidth=nK+5;
	subplotsx=nK+1;
	subplotsy=nUtop+1;   
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

	ymincurs = 0;
	ymaxcurs = 0;
	ymingrip = 0;
	ymaxgrip = 0;
	yminsp = 0;
	ymaxsp = 0;
	for ii = 1:nUtop
		i = topunits(ii);
		model = models{i};
		stimfilt = model.k;
		const = model.dc;
		sphist = model.ihbas*model.ih;
		couples = model.ihbas*model.ih2;
		spfilters = [couples(:,topunits(1:(ii-1))), couples(:,topunits((ii+1):end))];
		ymincurs = min(ymincurs, min(min(stimfilt(:,1:3))));
		ymaxcurs = max(ymaxcurs, max(max(stimfilt(:,1:3))));
		ymingrip = min(ymingrip, min(min(stimfilt(:,4))));
		ymaxgrip = max(ymaxgrip, max(max(stimfilt(:,4))));
		yminsp = min(yminsp, min(min(spfilters)));
		ymaxsp = max(ymaxsp, max(max(spfilters)));
	end

	for ii = 1:nUtop
		i = topunits(ii);
		model = models{i};
		idx = 1;
		stimfilt = model.k;
		const = model.dc;
		sphist = model.ihbas*model.ih;
		couples = model.ihbas*model.ih2;
		spfilters = [couples(:,topunits(1:(ii-1))), sphist, couples(:,topunits((ii+1):end))];

		ax=axes('position',sub_pos{1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		name = processed.unitnames{i};
		filt = sphist;
		filt = flipud(filt);
		dt_filt = model.dt*RefreshRate;
		tt = (0:length(filt)-1)*dt_filt;
		tt = (tt-max(tt))/RefreshRate/model.dt;
		plot(tt, filt, 'k');
		if length(tt) > 1
			xlim([min(tt) max(tt)]);
		end
		if i == nM
			title(name);
		end
		if i == 1
			xlabel('time (ms)');
		end

		for j = 1:(nUtop)
			ax=axes('position',sub_pos{j+1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			if j == ii
				plot(0,0)
				axis off
			else
				name = processed.unitnames{topunits(j)};
				filt = spfilters(:,j);
				filt = flipud(filt);
				dt_filt = model.dt*RefreshRate;
				tt = (0:length(filt)-1)*dt_filt;
				tt = (tt-max(tt))/RefreshRate/model.dt;
				plot(tt, filt);
				ylim([yminsp, ymaxsp])
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
		end

		for j = nUtop+(1:4)
			ax=axes('position',sub_pos{j+1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			name = names{j-nUtop+1};
			filt = stimfilt(:,j-nUtop);
			dt_filt = model.dt*RefreshRate;
			tt = (0:length(filt)-1)*dt_filt;
			if j == nUtop+4
				plot(tt, filt, 'r');
				ylim([ymingrip, ymaxgrip])		
			else
				plot(tt, filt, 'Color', [0 .5 0]);
				ylim([ymincurs, ymaxcurs])		
			end
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
		ax=axes('position',sub_pos{nK+1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
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