function jackknifecoherence_all(wd, fn_ins, fn_out)
	%Set wd to working directory
	%wd = './results_L1_stampede/';
	%wd = '/Users/dk/Dropbox/AFK_Neuron/Code_GLM_ben/monkeyresults/coherence_networkglm/';
	h = figure;
	g = figure;
	j = figure;
	cmap = jet(length(fn_ins));
	for iii = 1:length(fn_ins)
		fn_in = fn_ins{iii};
		load([wd fn_in])
		%proc_withheld 
		%For each unit
		coh_05 = zeros(nU,1);
		coh_cpl_05 = zeros(nU,1);
		coh_SE_05 = zeros(nU,1);
		coh_SE_cpl_05 = zeros(nU,1);
	
		maxt = 20000;
		for icell = 1:nU
		    truesp = {};
		    simsp = {};
		    simsp_cpl = {};
		    runs = size(proc_withheld.trialstartend,1);
		    K = [];
		    N_sample_max = 0;
		    uncoupled = load(['./results/GLM_cell_simulation_' num2str(goodunits(icell)) '.mat']);
		    nspks = 0;
		    for i = 1:runs
		        ts = proc_withheld.trialstartend(i, 1);
		        te = proc_withheld.trialstartend(i, 2);
		        if ts > maxt | te > maxt
		        	break
		        end
		        truesp{i} = proc_withheld.spiketrain(1:maxt,icell);
		        truesp{i} = truesp{i}(ts:te);
		        simsp_cpl{i} = Rt_glm{icell}(:);
		        simsp_cpl{i} = simsp_cpl{i}(ts:te);
		        simsp{i} = uncoupled.Rt_glm(:);
		        simsp{i} = simsp{i}(ts:te);
		        N_sample_max = max(te-ts+1, N_sample_max);
		        nspks = nspks + sum(truesp{i});
		    end
		    runs = i-1;
		    Pad = 2^(1+nextpow2(N_sample_max));  
		    f=linspace(0, 1, Pad)*RefreshRate;
	
		    Coh_cpl = (zeros(runs, Pad));
		    Coh = (zeros(runs, Pad));
		    Power = (zeros(runs, Pad));
		    Power_cpl = (zeros(runs, Pad));
		    Power_true = (zeros(runs, Pad));
		
		    for i = 1:runs
		        %Compute coherence of all runs
		        [coh_cpl, pwr_true, pwr_cpl, k] = coherence(truesp{i}, simsp_cpl{i}, N_sample_max);
    		    [coh, pwr1, pwr]  = coherence(truesp{i}, simsp{i}, N_sample_max);
    	    	Coh_cpl(i,:) = coh_cpl;
		        Coh(i,:) = coh;
    		    Power(i,:) = pwr;
    	    	Power_cpl(i,:) = pwr_cpl;
    	    	Power_true(i,:) = pwr_true;
		        K(i) = k;
		    end
		
			Power = mean(Power);
			Power_cpl = mean(Power_cpl);
			Power_true = mean(Power_true);
	
		    %Only compute for not NaN results
		    Coh = Coh(~isnan(Coh(:,1)),:);
		    if size(Coh,1) > 1
		        [Coh, Coh_SE] = jackknife_coh(Coh);
		    end
		
		    Coh_cpl = Coh_cpl(~isnan(Coh_cpl(:,1)),:);
		    if size(Coh_cpl,1) > 1
		        [Coh_cpl, Coh_cpl_SE] = jackknife_coh(Coh_cpl);
		    end
		    K = mean(K(~isnan(Coh(:,1))));
		    DF = runs*K-1; Null = sqrt(1-0.05^(1/DF));
		    L = 5; 
		    ii = 1:fix(Pad/L);
		    ii_s = 6;
	
		    coh1hz_cpl = abs(Coh_cpl);
		    coh1hz = abs(Coh);
		    meancoh1hz_cpl(icell) = mean(coh1hz_cpl(2:11));
		    meancoh1hz(icell) = mean(coh1hz(2:11));
	
			coh_05(icell) = Coh(ii_s);
			coh_cpl_05(icell) = Coh_cpl(ii_s);
			coh_SE_05(icell) = Coh_SE(ii_s);
			coh_SE_cpl_05(icell) = Coh_cpl_SE(ii_s);
		end

		%Plot scatter plot:
		figure(g);
		if iii == 1
			xlim([0 1])
			ylim([0 1])
			plot([0 1], [0 1], 'k')
			xlabel('Coherence b/w true spikes and uncoupled GLM')
			ylabel('Coherence b/w true spikes and coupled GLM')
		end
		hold on
		errorbar(coh_05, coh_cpl_05, coh_SE_cpl_05, '.', 'Color', cmap(iii, :));
		he = herrorbar(coh_05, coh_cpl_05, coh_SE_05, '.');
		set(he, 'MarkerEdgeColor', cmap(iii,:));

		%Plot scatter plot mean zero to one Hz coherence:
		figure(h);
		if iii == 1
			xlim([0 1])
			ylim([0 1])
			plot([0 1], [0 1], 'k')
			xlabel('Mean (0-1hz) Coherence b/w true spikes and uncoupled GLM')
			ylabel('Mean (0-1hz) Coherence b/w true spikes and coupled GLM')
		end
		hold on
		scatter(meancoh1hz, meancoh1hz_cpl, [], cmap(iii,:));
		plot(meancoh1hz, meancoh1hz_cpl, 'Color', cmap(iii,:));

		figure(j);
		if iii == 1
			%logl = logl_glm(3:end,:);
			cmap = jet(6);
			plot([0,10],[0 0], '--k')
			legend('\lambda = 1 (fully coupled)', '\lambda = 3', '\lambda = 10', '\lambda = 30', '\lambda = 100', '\lambda = 300', 'location', 'eastoutside')
			xlabel('Unit');
			ylabel('Coupled - uncoupled <|coh|> @ 1Hz')
		end
		hold all
	    plot(1:length(meancoh1hz), meancoh1hz_cpl - meancoh1hz, 'o', 'Color', cmap(iii,:))
	end
	saveplot(g, [wd fn_out 'coherence_scatter_all.eps'])
	saveplot(h, [wd fn_out '_coherence_scatter_mean1hz_all.eps'])
	saveplot(j, [wd fn_out '_coherence_scatter_mean1hz_all_relative.eps'], 'eps', [7 4])
end

function [Coh_all, Coh_SE] = jackknife_coh(Coh)
	nS = size(Coh,1);
	Coh_all = abs(mean(Coh));
	Coh_t = repmat(Coh_all, nS, 1);
	Coh_dropped = [];
	for i = 1:nS
		ii = [1:(i-1), (i+1):nS];
		Coh_dropped(i,:) = abs(mean(Coh(ii, :)));
	end
	Coh_SE = sqrt((nS-1)/nS*sum((Coh_dropped - Coh_t).^2,1));
end