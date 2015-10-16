function [processed, processed_withheld] = preprocess_crossval(datafile, binsize, dt, frames, folds, fold, stddev, goodunits)
	%Preprocess both spike data and stim data. Split data into 'folds' for cross-validation
	%
	%Usage:
	%		[processed, processed_withheld] = preprocess_crossval(datafile, binsize, dt, frames)
	%
	%Input:
	%		datafile = .mat file with data
	%		binsize = size in seconds of each time bin for stim
	%		dt = relative size of bins spike times are provided at		
	%		frames = number of previous stim frames to include
	%		folds = number of folds
	%		fold = fold indx
	%		stddev = (optional, default = 0) whether to normalize stim
	%		goodunits = (optional) if provided, will only return unit indices specified
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
	%	datafile = './mabel_reaching_5-4-10.mat';
	%	binsize = 1/100;
	%	dt = 0.1;
	%	frames = 4;
	%	folds = 10;
	%	fold = 3;
	%	processed = preprocess_crossval(datafile, binsize, dt, frames, folds, fold);

	%Event legend
	if (nargin < 5) stddev = 0; end
	TRIALEND = 15;
	TARGETAPPEARED = 110;
	TARGETREACHED = 111;
	ds = 0.01;

    load(datafile);
    processed.binsize = binsize;
    processed.cursor = [Cursor_X, Cursor_Y, Cursor_Z];
    processed.grip = Grip_force;
    processed.stim = [processed.cursor, processed.grip];
    %Normalize
    if stddev == 1
    	for idx = 1:size(processed.stim,2)
    		idx
	    	processed.stim(:,idx) = (processed.stim(:,idx)-mean(processed.stim(:,idx)))/std(processed.stim(:,idx));
    	end
	end
    processed.stim = resample(processed.stim, 1000*ds, 1000*binsize);
    processed.dt = dt;
    unitnames = who('CSPIK*');
    processed.unitnames = unitnames;
    processed.frames = frames;
    nU = length(unitnames);
    nB = size(processed.stim,1);
    nS = size(processed.stim,2);

    processed.spikes = {};
    processed.spiketrain = zeros(nB, nU);
    for idx = 1:nU
    	%Get data for idxth spike train
    	eval(['spikes = ' unitnames{idx} ';']);
    	processed.spikes{idx} = spikes*dt;
    	bins = ceil(spikes*dt);
    	for b = bins
	    	processed.spiketrain(b,idx) = processed.spiketrain(b,idx)+1;
		end
	end

    %Stack stim to include past and future frames relative to spike time
    processed.stacked = zeros(nB, nS*(2*frames+1));
    offsets = -frames:frames;
	nx = nS;
	nt = length(offsets);
	stim = vertcat(zeros(frames, nx), processed.stim, zeros(frames, nx));
	nB = size(stim,1);
	for i = 1:nx
	    for j = 1:nt
	        offset = offsets(j);
	        jj = (frames+offset+1):(nB-frames+offset);
	        processed.stacked(:,(i)*nt-j+1) = stim(jj, i);
	    end
	end

	%Here trial start is determined by Est_Movement_Initiation variable...
	%Note which bins are inside a trial
	trialendidx = find(Events_Data(2,:)==TRIALEND);
	trialstartbins = ceil(Est_Movement_Initiation*dt);
	trialendbins = ceil(Events_Data(1,trialendidx)*dt);
	%Only take the trial ends that are immediately after a movement start
	tendbins = [];
	for idx = 1:length(trialstartbins)
		startbin = trialstartbins(idx);
		%Find next trial end
		tend_idx = find(trialendbins > startbin, 1);
		tendbin = trialendbins(tend_idx);
		if ~isempty(tendbin)
			%If this is before the _next_ movement initiation, then take the pair
			nextstart_idx = find(trialstartbins > startbin,1);
			nextstart = trialstartbins(nextstart_idx);
			if isempty(nextstart)
				tendbins = [tendbins, tendbin];
			else 
				if nextstart > tendbin
					tendbins = [tendbins, tendbin];
				end
			end
		end
	end
	trialendbins = tendbins;

	trialchanges = zeros(size(processed.grip));
	trialchanges(trialstartbins) = 1;
	trialchanges(trialendbins) = -1;
	processed.intrial = cumsum(trialchanges);
	%Start outside of of a trial
	processed.intrial(1) = 0;
	processed.trialstartend = [trialstartbins', trialendbins'];

	%Split into folds here
	processed_withheld = processed;
	assert(fold <= folds, 'Must set test data fold to be less than or equal to fold number');
	blocksize = floor(nB/folds);
	preidx = 1:((fold-1)*blocksize);
	testidx = ((fold-1)*blocksize+1):(fold*blocksize);
	postidx = ((fold)*blocksize+1):nB;
	trainingidx = [preidx, postidx];
	trainmax = max(preidx);
	testmax = max(testidx);

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

	%Move test spikes relative to start of test data
	for idx = 1:nU
	    sp = processed.spikes{idx};
    	sptrainpre = sp(sp<trainmax);
    	sptrainpost = sp(sp>testmax);
	    sptest = sp(sp>trainmax);
    	processed.spikes{idx} = [sptrainpre, sptrainpost-testmax];
	    processed_withheld.spikes{idx} = sptest-trainmax;
	end

	processed_withheld.trialstartend = [];
	for idx = 1:size(processed.trialstartend,1)
		if processed.trialstartend(idx,1) > trainmax
			tstart = processed.trialstartend(idx,1)-trainmax;
			tend = processed.trialstartend(idx,2)-trainmax;
			processed_withheld.trialstartend = [processed_withheld.trialstartend; tstart, tend];
		elseif processed.trialstartend(idx,1) > testmax
			tstart = processed.trialstartend(idx,1)-testmax;
			tend = processed.trialstartend(idx,2)-testmax;
			processed.trialstartend(idx,:) = [tstart, tend];			
		end
	end

	if nargin == 6
		[processed, processed_withheld] = remove_units(goodunits, processed, processed_withheld);
	end
end


