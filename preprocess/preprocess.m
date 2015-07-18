function processed = preprocess(stimfile, responsefiles, binsize, unitidx)
	%Preprocess both spike data and stim data. To save memory for later fitting,
	%only load the stimulus of one unit at a time. 
	%
	%Usage:
	%		processed = preprocess(stimfile, responsefile, binsize, unitidx)
	%
	%Input:
	%		stimfile = whitenoise input data
	%		responsefiles = list of isk files containing response spike trains
	%		binsize = size in seconds of each time bin
	%		unitidx = index of element in response files whose stim to load
	%	
	%Output:
	%		processed is a structure containing the following fields:
	%			spikes = [nB x nU] array with spikes from all channels binned 
	%				according to binsize. 
	%				nB = no. bins, nU = no. units.
	%			stim = [nB x (nX x nX)] array of stimuli presented to 'unitidx'th cell
	%			binsize = binsize used
	%			unitnames = names of units loaded into spikes
	%			
	%
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1/30;
	%	unitidx = 2;
	%	processed = preprocess(stimfile, rf, binsize, unitidx);
	
	%Long params
	iL = 2;
	stim_length = {'short', 'long'};
    load(['RetinaCellParameters_' stim_length{iL} '.mat']) ;

    processed.binsize = binsize;
    nU = length(responsefiles);
    maxlag = max(lagshifts);
    unitnames = {};
    for idx = 1:length(responsefiles)
    	if isstruct(responsefiles)
	    	responsefile = responsefiles(idx).name;
	    elseif iscell(responsefiles)
	    	responsefile = responsefiles{idx};
	    end
    	icell = str2num(responsefile(end-5:end-4));
    	if isempty(icell)
    		icell = str2num(responsefile(end-4));
    	end
    	unitnames = {unitnames{:}, num2str(icell)};
	    %Responses
	    lag = lagshifts(icell) ;
    	fsize = Nv(icell)^2;
    	fres = fopen(responsefile, 'rt');
    	R = textscan(fres, '%u\n');
    	fclose(fres);
    	R = double(R{1,1});
    	T = length(R);
    	if ~isfield(processed, 'spikes')
    		processed.spikes = zeros(T, nU);
    	end
	    
	    %Stimulus
	    if idx == unitidx
	    	fstim = fopen(stimfile, 'rb');
	    	S = ReadFramev2(fstim,T,Nx,Nv(icell),cx,x0(icell),y0(icell));
	    	fclose(fstim);	
	    	pX = length(S)/T;
    		S = reshape(S, [pX T])';
	    	%Remove mean, standardize by std 
    		S = S - repmat(mean(S), [T 1]);
    		S = S./repmat(std(S), [T 1]);
	    	%Save in processed structure
    		processed.stim = S;
    		processed.unit = num2str(icell);
    		processed.Nv = Nv(icell);
	    end
   	    %Add lag to some cells
	    R = circshift(R,-lag);
    	%Save in processed structure
    	processed.spikes(:,idx) = R;
	end

	%Truncate to max lag
    processed.stim = processed.stim(1:end-maxlag,:,:);
    processed.spikes = processed.spikes(1:end-maxlag,:);
    processed.unitidx = unitidx;
    processed.unitnames = unitnames;

