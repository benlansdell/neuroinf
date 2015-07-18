function data = filters_monkey_sp_stm(processed, nK_sp, nK_stm, dt_sp, dt_stm)
	%Prepare spike and stimulus data for GLM
	%
	%Usage:
	%	data = filters_monkey_sp_stm(processed, nK_sp, nK_stm)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	nK_sp = number of timebins used for spike history filter
	%	nK_stm = number of timebins used for stimulus filter
	%	dt_sp = (optional, default = binsize in processed structure) step size of spike history filter
	%		in seconds. Must be a multiple of the data's binsize.
	%	dt_pos = (optional, default = binsize in processed structure) step size of position filter in
	%		seconds. Must be a multiple of the data's binsize	
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [1 x nB] array where y_i is the number of spikes at time bin i
	%		X = [nB x nK] array where X_ij is the value of covariate j, at
	%			time bin i. Note: nK = nK_sp + nK_stm
	%		k = Names of each filter, a [n x 3] cell array in which each row is 
	%			of the form
	%				['filter j name', [idxj1 idxj2 ...], binsize] 
	%			where the second column lists indices in 1:nK which belong to the
	%			filter, the third column indicates the temporal resultion that filter
	%			runs at
	%		stim = stim data trimmed in the same way X and y are. 
	%		sp_hist = indices corresponding to (auto) spike history terms
	%		nK_sp = number of spike history components
	%		nK_stm = number of stimulus components
	%
	%Test code:
	%	datafile = './data/mabel_reaching_5-4-10.mat';
	%	binsize = 1/100;
	%	nK_sp = 20;
	%	nK_stm = 6;
	%	dt_sp = binsize;
	%	dt_stm = 5/100;
	%	unitidx = 13;
	%	processed = preprocess_monkey(datafile, binsize, unitidx);
	%	data = filters_monkey_sp_stm(processed, nK_sp, nK_stm, dt_sp, dt_stm);

	if (nargin < 4) dt_sp = processed.binsize; end
	if (nargin < 5) dt_pos = processed.binsize; end

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_stm,processed.binsize)==0, 'Invalid dt_stm. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_stm = dt_stm/processed.binsize;

	nB = size(processed.spikes,1);
	nU = 1;
	unit = processed.unitidx;
	nK = nU*nK_sp + 4*nK_stm;

	data.X = zeros(nB, nK);
	data.k = cell(5,3);
	data.k{1,1} = 'spike history'; 
	data.k{1,2} = 1:nK_sp;
	data.k{1,3} = dt_sp;
	data.k{2,1} = 'curs x';
	data.k{2,2} = (1:(nK_stm)) + nK_sp;
	data.k{2,3} = dt_stm;
	data.k{3,1} = 'curs y';
	data.k{3,2} = (1:(nK_stm)) + nK_sp + nK_stm;
	data.k{3,3} = dt_stm;
	data.k{4,1} = 'curs z';
	data.k{4,2} = (1:(nK_stm)) + nK_sp + 2*nK_stm;
	data.k{4,3} = dt_stm;
	data.k{5,1} = 'grip';
	data.k{5,2} = (1:(nK_stm)) + nK_sp + 3*nK_stm;
	data.k{5,3} = dt_stm;
	%Record specifically which indices are spike history indices for model simulation
	data.sp_hist = data.k{1,2};

	%Make stimulus vector at each timebin
	for j = (nK_sp*steps_sp+1:nB-nK_stm*steps_stm)
		%(past) spike history
		shist = processed.spikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, unit);
		%(future) cursor
		curx = processed.cursor(j:steps_stm:(j+(nK_stm-1)*steps_stm),1);
		cury = processed.cursor(j:steps_stm:(j+(nK_stm-1)*steps_stm),2);
		curz = processed.cursor(j:steps_stm:(j+(nK_stm-1)*steps_stm),3);
		%(future) grip
		grip = processed.grip(j:steps_stm:(j+(nK_stm-1)*steps_stm),1);
		%Form stim vector
		data.X(j,:) = [shist' curx' cury' curz' grip'];
	end

	%Truncate to exclude start of recording where spike history isn't well defined
	data.X = data.X((nK_sp*steps_sp+1):(nB-nK_stm*steps_stm),:);
	data.y = processed.spikes((nK_sp*steps_sp+1):(nB-nK_stm*steps_stm), unit)';
	data.cursor = processed.cursor((nK_sp*steps_sp+1):(nB-nK_stm*steps_stm),:); 
	data.grip = processed.grip((nK_sp*steps_sp+1):(nB-nK_stm*steps_stm),:);

	data.nK_sp = nK_sp;
	data.nK_stm = nK_stm;