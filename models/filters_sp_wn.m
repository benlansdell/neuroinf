function data = filters_sp_stm(processed, nK_sp, nK_stm, dt_sp, dt_stm)
	%Prepare spike and stimulus data for GLM
	%includes spike history and stimulus filters:
	%
	%	y(i) ~ Pn(g(eta_i))
	%
	%where
	%
	%	eta_i = \sum y(i-j) k_y(i) + \sum x(i-j) k_x(j)
	%
	%Usage:
	%	data = filters_sp_stm(processed, nK_sp, nK_stm, dt_sp, dt_stm)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	nK_sp = number of timebins used for spike history filter
	%	nK_stm = number of timebins used for stimulus filter
	%	dt_sp = (optional, default = binsize in processed structure) step size of spike history filter
	%		in seconds. Must be a multiple of the data's binsize.
	%	dt_stm = (optional, default = binsize in processed structure) step size of position filter in
	%		seconds. Must be a multiple of the data's binsize
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [nU x nB] array where y_ij is the number of spikes at time bin j for unit i.
	%		X = [nU x nB x nK] array where X_ijk is the value of covariate k, at time bin j, for unit i
	%			Note: nK = nK_sp + 2*nK_pos
	%		k = Names of each filter, a [n x 2] cell array in which each row is of the form ['filter j name', [idxj1 idxj2 ...]]
	%			Note: The second column lists indices in 1:nK to which the label applies
	%		torque = torque data trimmed in the same way X and y are. 
	%			Note: truncated at start and end because spike and cursor trajectory are not defined for first 
	%			and last nK_sp and nK_pos timebins respectively.
	%		dtorque = trimmed dtorque
	%		ddtorque = trimeed ddtorque
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);

	if (nargin < 4) dt_sp = processed.binsize; end
	if (nargin < 5) dt_pos = processed.binsize; end

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_pos,processed.binsize)==0, 'Invalid dt_pos. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_pos = dt_pos/processed.binsize;

	nU = size(processed.binnedspikes,2);
	nB = size(processed.binnedspikes,1);
	nK = nK_sp + 2*nK_pos;

	data.X = zeros(nU, nB, nK);
	data.k = cell(3,3);
	data.k{1,1} = 'spike history'; 
	data.k{1,2} = 1:nK_sp;
	data.k{1,3} = dt_sp;
	data.k{2,1} = 'RU pos'; 
	data.k{2,2} = (1:nK_pos) + nK_sp;
	data.k{2,3} = dt_pos;
	data.k{3,1} = 'FE pos'; 
	data.k{3,2} = (1:nK_pos) + nK_sp + nK_pos;
	data.k{3,3} = dt_pos;
	%Record specifically which indices are spike history indices for model simulation
	data.sp_hist = data.k{1,2};

	%For each unit, add data to X array
	for idx=1:nU 
		%Make stimulus vector at each timebin
		for j = (nK_sp*steps_sp+1):(nB-nK_pos*steps_pos)
			%(past) spike history
			shist = processed.binnedspikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx);
			%(future) torque trajectory
			torqueRU = processed.torque(j:steps_pos:(j+(nK_pos-1)*steps_pos),1);
			torqueFE = processed.torque(j:steps_pos:(j+(nK_pos-1)*steps_pos),2);
			%Add a small amount of normal noise to torque data to prevent rank deficient matrices...
			%torqueRU = torqueRU + randn(size(torqueRU))/10;
			%torqueFE = torqueFE + randn(size(torqueFE))/10;
			%Form stim vector
			data.X(idx,j,:) = [shist' stim'];
		end
	end
	%Truncate to exclude start and end of recording where spike history 
	%and cursor trajectory aren't well defined
	data.X = data.X(:,(nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:); %(nkt+1:end-nkt,:);
	data.y = processed.binnedspikes((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos), :)';
	%Truncate other data for comparison, too
	data.torque = processed.torque((nK_sp*steps_sp+1):(nB-nK_pos*steps_pos),:); 
