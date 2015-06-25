function data = filters_sprc_stm(processed, nK_sp, nK_stm)
	%Prepare spike and stimulus data for GLM
	%includes spike history (in raised cosine basis) and stimulus filters:
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
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [nU x nB] array where y_ij is the number of spikes at time bin j
	%			for unit i.
	%		X = [nU x nB x nK] array where X_ijk is the value of covariate k, at
	%			time bin j, for unit i. Note: nK = nK_sp + nK_stm
	%		k = Names of each filter, a [n x 2] cell array in which each row is 
	%			of the form ['filter j name', [idxj1 idxj2 ...]] 
	%			Note: The second column lists indices in 1:nK to which the label
	%			applies
	%		stim = stim data trimmed in the same way X and y are. 
	%
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 6;
	%	nK_stm = 6;
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm(processed, nK_sp, nK_stm);


	dt_sp = processed.binsize;
	dt_stm = processed.binsize;
	unit = processed.unitidx;

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_stm,processed.binsize)==0, 'Invalid dt_stm. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_stm = dt_stm/processed.binsize;

	nB = size(processed.stim,1);
	nU = 1;
	nS = size(processed.stim,2);

	T = nK_sp*dt_sp;
	[rcbasis, spbasis, nK_rc] = makeRCBasis(dt_sp, T);
	nK = nK_rc + nS*nK_stm;

	data.X = zeros(nB, nK);
	data.k = cell(2,3);
	data.k{1,1} = 'spike history'; 
	data.k{1,2} = 1:nK_rc;
	data.k{1,3} = dt_sp;
	data.k{2,1} = 'stim';
	data.k{2,2} = (1:(nK_stm*nS)) + nK_rc;
	data.k{2,3} = dt_stm;
	%Record specifically which indices are spike history indices for model simulation
	data.sp_hist = data.k{1,2};

	%Make stimulus vector at each timebin
	for j = (nK_sp*steps_sp+1:nB)
		%(past) spike history
		%project_rc(processed.binnedspikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx), rcbasis);
		shist = project_rc(processed.spikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, unit), rcbasis);
		%(past) stimulus
		stim = processed.stim(j-nK_stm*steps_stm:steps_stm:j-steps_stm,:);
		stim = reshape(squeeze(stim)', 1, []);
		%Form stim vector
		data.X(j,:) = [shist' stim];
	end
	%Truncate to exclude start of recording where spike history isn't well defined
	data.X = data.X((nK_sp*steps_sp+1):end,:); %(nkt+1:end-nkt,:);
	data.y = processed.spikes((nK_sp*steps_sp+1):end, unit)';
	%Truncate other data for comparison, too
	data.stim = processed.stim((nK_sp*steps_sp+1):end,:); 
	data.nK_sp = nK_sp; 
	data.nK_stm = nK_stm;

	data.rcbasis = rcbasis;
	data.spbasis = spbasis;
end

function sphistory_rc = project_rc(sphistory, rcbasis)
	sphistory_rc = rcbasis*sphistory;
end

function [rcbasis, spbasis, nK_rc] = makeRCBasis(dt, T)
	%Define basis of raised cosine functions
	%Hard-wired log-scale
	a = 15;
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