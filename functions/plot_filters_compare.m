function plot_filters_compare(allmodels, processed, fn_out, refreshrates)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	models = data structure output by function in ./fitting (containing fitted coefficients)
	%	fn_out = base filename to write plots to for each unit

	%Plot colors:

	if nargin < 4
		refreshrates = 100*ones(length(allmodels), 1);
	end
	cmap = get(groot, 'DefaultAxesColorOrder');

	nG = length(allmodels);
	global RefreshRate;
	nU = 1; %number of units
	nM = length(allmodels{1});
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

	%Before plotting, collect the largest cursor and grip values to set scale by
	ymincurs = 0;
	ymaxcurs = 0;
	ymingrip = 0;
	ymaxgrip = 0;
	for j = 1:nG
		models = allmodels{j};
		for i = 1:nM
			model = models{i};
			idx = 1;
			stimfilt = model.k;
			stimfilt = reshape(stimfilt, [], 4);
			%if 
			curs = stimfilt(:,1:3);
			grip = stimfilt(:,4);
			ymincurs = min(ymincurs, min(min(curs)));
			ymaxcurs = max(ymaxcurs, max(max(curs)));
			ymingrip = min(ymingrip, min(grip));
			ymaxgrip = max(ymaxgrip, max(grip));
		end
	end
	clear stimfilt;

	clf
	for i = 1:nM
		for j = 1:nG
			models = allmodels{j};
			model = models{i};
			idx = 1;
			stimfilt{j} = model.k;
			stimfilt{j} = reshape(stimfilt{j}, [], 4);
			const = model.dc;
			sphist{j} = model.ihbas*model.ih;
		end

		for j = 1:(nK)
			%Extract data
			name = names{j};

			if j == 1
				for g = 1:nG
					filt{g} = sphist{j};
					filt{g} = flipud(filt{g});
				end
				%ymin = min([min(filt), min(filt2)])*1.2;
				%ymax = max([max(filt), max(filt2)])*1.2;
			elseif j > 1 & j < 5
				ymin = ymincurs*1.2;
				ymax = ymaxcurs*1.2;
				for g = 1:nG
					filt{g} = stimfilt{g}(:,j-1);
				end
			else 
				ymin = ymingrip*1.2;
				ymax = ymaxgrip*1.2;
				for g = 1:nG
					filt{g} = stimfilt{g}(:,j-1);
				end
			end
			ax=axes('position',sub_pos{j,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			if j > 1
				ylim([ymin ymax])
			end
			hold on
			for g = 1:nG
				clr = cmap(g,:);
				models = allmodels{g};
				RefreshRate = refreshrates(g);
				model = models{i};
				dt_filt = model.dt*RefreshRate;
				tt = (0:length(filt{g})-1);%*dt_filt;
				if j == 1
					tt = (tt-max(tt));
				else
					tt = (tt-max(tt)/2)*1000/RefreshRate;
				end
				plot(tt, filt{g}, 'Color', clr);
			end
			%if length(tt) > 1
			%	xlim([min(tt) max(tt)]);
			%end
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
		str1(4) = {['Number of spikes: ' num2str(sum(processed.spiketrain(:,i)))]};
		text(0.1,0.8,str1, 'FontSize', 5, 'Interpreter', 'none')
		axis off
	end
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);	
end