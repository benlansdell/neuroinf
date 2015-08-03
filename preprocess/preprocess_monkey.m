function processed = preprocess_monkey(datafile, binsize, unitidx)
	%Preprocess both spike data and stim data
	%
	%Usage:
	%		processed = preprocess_monkey(datafile, binsize)
	%
	%Input:
	%		datafile = .mat file with data
	%		binsize = size in seconds of each time bin
	%		unitidx = index of element in response files whose stim to load
	%	
	%Output:
	%		processed is a structure containing the following fields:
	%			spikes = [nB x nU] array with spikes from all channels binned 
	%				according to binsize. 
	%				nB = no. bins, nU = no. units.
	%			cursor = [nB x 3] array of stimuli presented to 'unitidx'th cell
	%			grip = [nB x 1] array of grip force
	%			binsize = binsize used
	%			unitnames = names of units loaded into spikes
	%			unitidx = index of unit to fit GLM to
	%
	%Test code:
	%	datafile = './data/mabel_reaching_5-4-10.mat';
	%	binsize = 1/100;
	%	unitidx = 13;
	%	processed = preprocess_monkey(datafile, binsize, unitidx);

	%Event legend
	TRIALSTART = 10;
	TRIALEND = 15;

    load(datafile) ;
    processed.binsize = binsize;
    processed.cursor = [Cursor_X, Cursor_Y, Cursor_Z];
    processed.grip = Grip_force;
    unitnames = who('CSPIK*');
    processed.unitnames = unitnames;
    processed.unitidx = unitidx;
    nU = length(unitnames);
    nB = size(processed.cursor,1);
    %nB = ceil(nB/3);
    processed.spikes = zeros(nB, nU);
    for idx = 1:nU
    	%Get data for idxth spike train
    	eval(['spikes = ' unitnames{idx} ';']);
    	%Bin in processed.spikes matrix @ 100Hz
    	bins = ceil(spikes/10);
    	for b = bins
	    	processed.spikes(b, idx) = processed.spikes(b, idx) + 1;
		end
	end

	%Note which bins are inside a trial
	trialstartidx = find(Events_Data(2,:)==TRIALSTART);
	trialendidx = find(Events_Data(2,:)==TRIALEND);
	trialstartbins = ceil(Events_Data(1,trialstartidx)/10);
	trialendbins = ceil(Events_Data(1,trialendidx)/10);
	trialchanges = zeros(size(processed.grip));
	trialchanges(trialstartbins) = 1;
	trialchanges(trialendbins) = -1;
	processed.intrial = cumsum(trialchanges);

	%Remove start of recording, since cursor and grip are zero
	startidx = find(processed.grip>0,1);
	processed.cursor = processed.cursor(startidx:end,:);
	processed.grip = processed.grip(startidx:end,:);
	processed.spikes = processed.spikes(startidx:end,:);
	processed.intrial = processed.intrial(startidx:end,:);
	%Start outside of of a trial
	processed.intrial(1) = 0;
	processed.trialstartend = [trialstartbins'-startidx+1, trialendbins'-startidx+1];