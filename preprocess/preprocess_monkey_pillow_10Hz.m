function processed = preprocess_monkey_pillow_10Hz(datafile, binsize, dt, frames)
	%Preprocess both spike data and stim data
	%
	%Usage:
	%		processed = preprocess_monkey_pillow(datafile, binsize, dt, unitidx)
	%
	%Input:
	%		datafile = .mat file with data
	%		binsize = size in seconds of each time bin for stim
	%		dt = relative size of bins spike times are provided at		
	%		frames = number of previous stim frames to include
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
	%			stacked = stacked stim for computing STA
	%
	%Test code:
	%	datafile = './data/mabel_reaching_5-4-10.mat';
	%	binsize = 1/10;
	%	dt = 0.01;
	%	frames = 4;
	%	processed = preprocess_monkey_pillow(datafile, binsize, dt, frames);

	%Event legend
	TRIALSTART = 10;
	TRIALEND = 15;

    load(datafile) ;
    processed.binsize = binsize;
    processed.cursor = [Cursor_X, Cursor_Y, Cursor_Z];
    processed.grip = Grip_force;
    processed.stim = [processed.cursor, processed.grip];
    processed.dt = dt;
 
    %Downsample stim
    every = 100*binsize;
    nB = size(processed.cursor,1);
    bins = 1:every:nB;
    processed.stim = processed.stim(bins,:);
    processed.grip = processed.grip(bins,:);
    processed.cursor = processed.cursor(bins,:);

    unitnames = who('CSPIK*');
    processed.unitnames = unitnames;
    processed.frames = frames;
    nU = length(unitnames);
    nB = size(processed.cursor,1);
    nS = size(processed.stim,2);

    processed.spikes = {};
    processed.spiketrain = zeros(nB, nU);
    for idx = 1:nU
    	%Get data for idxth spike train
    	eval(['spikes = ' unitnames{idx} ';']);
	    %Change scale spikes are specified at
    	processed.spikes{idx} = spikes*dt;
    	bins = ceil(spikes*dt);
    	for b = bins
	    	processed.spiketrain(b,idx) = processed.spiketrain(b,idx)+1;
		end
	end

    %Stack stim to include future frames from current one
    processed.stacked = zeros(nB, nS*frames)
	for idx = 1:frames
		processed.stacked(1:(end-frames), (idx-1)*nS+1:idx*nS) = processed.stim(1+idx-1:end-frames+idx-1,:);
	end
	%Trim end of stim, stacked, and offset spike times
	processed.stacked = processed.stacked(1:end-frames,:);
	processed.stim = processed.stim(1:end-frames,:);
	processed.cursor = processed.cursor(1:end-frames,:);
	processed.grip = processed.grip(1:end-frames,:);
	processed.spiketrain = processed.spiketrain(1:end-frames,:);

	%Note which bins are inside a trial
	trialstartidx = find(Events_Data(2,:)==TRIALSTART);
	trialendidx = find(Events_Data(2,:)==TRIALEND);
	trialstartbins = ceil(Events_Data(1,trialstartidx)*dt);
	trialendbins = ceil(Events_Data(1,trialendidx)*dt);
	trialchanges = zeros(size(processed.grip));
	trialchanges(trialstartbins) = 1;
	trialchanges(trialendbins) = -1;
	processed.intrial = cumsum(trialchanges);

	%Remove start of recording, since cursor and grip are zero
	startidx = find(processed.grip>0,1);
	processed.stacked = processed.stacked(startidx:end,:);
	processed.stim = processed.stim(startidx:end,:);
	processed.cursor = processed.cursor(startidx:end,:);
	processed.grip = processed.grip(startidx:end,:);
	processed.spiketrain = processed.spiketrain(startidx:end,:);
	for idx = 1:nU
		processed.spikes{idx} = processed.spikes{idx}-(startidx-1);
		processed.spikes{idx} = processed.spikes{idx}(processed.spikes{idx}>0);
	end
	processed.intrial = processed.intrial(startidx:end,:);
	%Start outside of of a trial
	processed.intrial(1) = 0;
	processed.trialstartend = [trialstartbins'-startidx+1, trialendbins'-startidx+1];