function [MSTM, SPNDS] = trimextratrial(MSTM, SPNDS, processed)
	nB = size(MSTM,1);
	pt = [];
	dt = processed.dt;
	for idx = 1:size(processed.trialstartend,1)
		tstart = processed.trialstartend(idx,1);
		tend = processed.trialstartend(idx,2);
		if tend < nB
			pt(idx,1:2) = [tstart, tend];
		else
			break
		end
	end
	%Chop out spikes past end of stim
	%sp = SPNDS;
	SPNDS(SPNDS*dt>nB) = [];
	%Make list of intervals to be excluded from data
	trialendstart = [[1; pt(:,2)], [pt(:,1); nB]];
	%Chop out spikes that are outside of a trial
	for idx = 1:190 %size(trialendstart,1)
		trimstart = trialendstart(idx,1);
		trimend = trialendstart(idx,2);
		MSTM(trimstart:trimend,:) = [];
		Dt = (trimend-trimstart+1)/processed.dt;
		SPNDS((SPNDS*dt>trimstart) & (SPNDS*dt<trimend)) = [];
		SPNDS(SPNDS*dt>trimend) = SPNDS(SPNDS*dt>trimend)-Dt;
	end
end

