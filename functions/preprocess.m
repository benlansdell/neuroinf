function [processed, processed_withheld] = preprocess(datafile, binsize, dt, frames)
	%Preprocess both spike data and stim data
	%
	%Usage:
	%		[processed, processed_withheld] = preprocess(datafile, binsize, dt, frames)
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
	%		processed_withheld is a structure containing 20% of data withheld for
	%			testing
	%
	%Test code:
	%	datafile = './data/mabel_reaching_5-4-10.mat';
	%	binsize = 1/100;
	%	dt = 0.1;
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
    	processed.spikes{idx} = spikes*dt;
    	bins = ceil(spikes/10);
    	for b = bins
	    	processed.spiketrain(b,idx) = processed.spiketrain(b,idx)+1;
		end
	end

    %Stack stim to include past and future frames relative to spike time
    processed.stacked = zeros(nB, nS*(2*frames+1));
    offsets = -frames:frames;
	for idx = 1:length(offsets)
		offset = offsets(idx);
		processed.stacked(1+frames:(end-frames), ((idx-1)*nS+1):idx*nS) = ...
		 processed.stim(1+frames+offset:end-frames+offset,:);
	end

	%Note which bins are inside a trial
	trialstartidx = find(Events_Data(2,:)==TRIALSTART);
	trialendidx = find(Events_Data(2,:)==TRIALEND);
	trialstartbins = ceil(Events_Data(1,trialstartidx)*dt);
	trialendbins = ceil(Events_Data(1,trialendidx)*dt);
	trialchanges = zeros(size(processed.grip));
	trialchanges(trialstartbins) = 1;
	trialchanges(trialendbins) = -1;
	processed.intrial = cumsum(trialchanges);
	%Start outside of of a trial
	processed.intrial(1) = 0;

	%Split into 20% test and 80% training set
	nB = size(processed.spiketrain,1);
	trainmax = ceil(nB*0.8);
	processed_withheld = processed;
	trainingidx = (1:nB)<=trainmax;
	processed.cursor = processed.cursor(trainingidx,:);
	processed.grip = processed.grip(trainingidx,:);
	processed.stim = processed.stim(trainingidx,:);
	processed.spiketrain = processed.spiketrain(trainingidx,:);
	processed.stacked = processed.stacked(trainingidx,:);
	processed.intrial = processed.intrial(trainingidx);

	processed_withheld.cursor = processed_withheld.cursor(~trainingidx,:);
	processed_withheld.grip = processed_withheld.grip(~trainingidx,:);
	processed_withheld.stim = processed_withheld.stim(~trainingidx,:);
	processed_withheld.spiketrain = processed_withheld.spiketrain(~trainingidx,:);
	processed_withheld.stacked = processed_withheld.stacked(~trainingidx,:);
	processed_withheld.intrial = processed_withheld.intrial(~trainingidx);

	for idx = 1:nU
	    sp = processed.spikes{idx};
    	sptrain = sp(sp<trainmax);
	    sptest = sp(sp>trainmax);
    	processed.spikes{idx} = sp;
	    processed_withheld.spikes{idx} = sptest-trainmax;
	end
