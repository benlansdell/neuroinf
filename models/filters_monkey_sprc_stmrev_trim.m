function data = filters_monkey_sprc_stm_trim(processed, nK_sp, nK_stm, a, dt_sp, dt_stm)
	%Prepare spike and stimulus data for GLM
	%
	%Usage:
	%	data = filters_monkey_sp_stm_trim(processed, nK_sp, nK_stm)
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
	%	a = 10;
	%	dt_sp = binsize;
	%	dt_stm = 5/100;
	%	unitidx = 13;
	%	processed = preprocess_monkey(datafile, binsize, unitidx);
	%	data = filters_monkey_sprc_stm_trim(processed, nK_sp, nK_stm, a, dt_sp, dt_stm);

	if (nargin < 4) a = 15; end
	if (nargin < 5) dt_sp = processed.binsize; end
	if (nargin < 6) dt_pos = processed.binsize; end

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_stm,processed.binsize)==0, 'Invalid dt_stm. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_stm = dt_stm/processed.binsize;

	T = nK_sp*dt_sp;
	[rcbasis, spbasis, nK_rc] = makeRCBasis(dt_sp, T, a);

	nB = size(processed.spikes,1);
	nU = 1;
	unit = processed.unitidx;
	nKs = 2*nK_stm-1;
	nK = nU*nK_rc + 4*nKs;

	data.X = zeros(nB, nK);
	data.k = cell(5,3);
	data.k{1,1} = 'spike history'; 
	data.k{1,2} = 1:nK_rc;
	data.k{1,3} = dt_sp;
	data.k{2,1} = 'curs x';
	data.k{2,2} = (1:(nKs)) + nK_rc;
	data.k{2,3} = dt_stm;
	data.k{3,1} = 'curs y';
	data.k{3,2} = (1:(nKs)) + nK_rc + nKs;
	data.k{3,3} = dt_stm;
	data.k{4,1} = 'curs z';
	data.k{4,2} = (1:(nKs)) + nK_rc + 2*nKs;
	data.k{4,3} = dt_stm;
	data.k{5,1} = 'grip';
	data.k{5,2} = (1:(nKs)) + nK_rc + 3*nKs;
	data.k{5,3} = dt_stm;
	%Record specifically which indices are spike history indices for model simulation
	data.sp_hist = data.k{1,2};

	%Make stimulus vector at each timebin
	strpt = 1+max(nK_sp*steps_sp,(nK_stm-1)*steps_stm);
	for j = (strpt:nB-(nK_stm-1)*steps_stm)
		%(past) spike history
		shist = project_rc(processed.spikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, unit), rcbasis);
		%(past and future) cursor
		curx = processed.cursor(j-(nK_stm-1)*steps_stm:steps_stm:(j+(nK_stm-1)*steps_stm),1);
		cury = processed.cursor(j-(nK_stm-1)*steps_stm:steps_stm:(j+(nK_stm-1)*steps_stm),2);
		curz = processed.cursor(j-(nK_stm-1)*steps_stm:steps_stm:(j+(nK_stm-1)*steps_stm),3);
		%(past and future) grip
		grip = processed.grip(j-(nK_stm-1)*steps_stm:steps_stm:(j+(nK_stm-1)*steps_stm),1);
		%Form stim vector
		data.X(j,:) = [shist' curx' cury' curz' grip'];
	end

	%Trim to only include the trial times...
	trimstart = 1+max(nK_sp*steps_sp,(nK_stm-1)*steps_stm);
	trimend = (nK_stm-1)*steps_stm;
	trimmedpairs = [processed.trialstartend(:,1)+trimstart, processed.trialstartend(:,2)-trimend];
	validtrials = trimmedpairs(:,1) < trimmedpairs(:,2);
	trimmedpairs = trimmedpairs(validtrials,:);
	intrialchange = zeros(size(data.X,1),1);
	intrialchange(trimmedpairs(:,1)) = 1;
	intrialchange(trimmedpairs(:,2)) = -1;
	intrialidx = cumsum(intrialchange);

	data.X = data.X(intrialidx==1,:);
	data.y = processed.spikes(intrialidx==1, unit)';
	data.cursor = processed.cursor(intrialidx==1,:); 
	data.grip = processed.grip(intrialidx==1,:);

	data.nK_sp = nK_sp; 
	data.nK_rc = nK_rc;
	data.nK_stm = nK_stm;

	data.rcbasis = rcbasis;
	data.spbasis = spbasis;
end

function sphistory_rc = project_rc(sphistory, rcbasis)
	sphistory_rc = rcbasis*sphistory;
end

function [rcbasis, spbasis, nK_rc] = makeRCBasis(dt, T, a)
	%Define basis of raised cosine functions
	nTotal = floor(T/dt);
	%Create basis function function
	B = @(t, a, psi, phi) iif(a*log(t-psi)>phi-pi & a*log(t-psi)<phi+pi, 0.5*(cos((a*log(t-psi)-phi))+1), ...
		true, zeros(size(t)));
	tt = 0:dt:(T-dt);
	nT = length(tt);
	%Compute phi0 and psi for desired basis
	phi1 = a*log(dt*(1-exp(-pi/a))^(-1));
	phi0 = phi1-pi;
	psi = -exp(phi0/a);
	%Compute each function for phi_i
	phi = phi0;
	%Compute matrix of basis vectors
	rc = zeros(nT, nTotal);
	for j = 1:nTotal
		for i = 1:nT
			rc(i, j) = B(tt(i), a, psi, phi);
		end
		%Normalize each column
		if norm(rc(:,j))>0
			rc(:,j) = rc(:,j)/norm(rc(:,j));
		end
		phi = phi + pi;
	end
	%Truncate to just the non-zero columns
	nK_rc = rank(rc);
	rc = rc(:,1:nK_rc);
	%Flip so time near spike is best resolved
	rc = flipud(rc);
	%Compute pseudo inverse of this matrix
	rcbasis = inv(rc'*rc)*rc';
	spbasis = rc;
end