function processed = preprocess_monkey_pillow(datafile, binsize, dt, frames)
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
	%	binsize = 1/100;
	%	dt = 0.1;
	%	frames = 4;
	%	processed = preprocess_monkey_pillow(datafile, binsize, dt, frames);

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
	%newend = ;
	%for idx = 1:nU
	%	processed.spikes{idx} = processed.spikes{idx}(processed.spikes{idx}<newend);
	%end

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
