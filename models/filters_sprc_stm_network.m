function data = filters_sprc_stm_network(processed, nK_sp, nK_stm, a)
	%Prepare spike and stim data for GLM which includes spike history and stim
	%history filters 
	%
	%	y(i) ~ Pn(g(eta_i))
	%
	%where
	%
	%	eta_i = \sum y(i-j) k_sp(i) + \sum x_1(i+j) k_1(j) + \sum x_2(i+j) k_2(j)
	%
	%Spike history filter is saved as a rasied cosine function coefficients
	%
	%Usage:
	%	data = filters_sprc_stm_network(processed, nK_sp, nK_stm)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	nK_sp = number of timebins used for spike history filter for all units
	%	nK_stm = number of timebins used for stim filter
	%	a = (optional, default = 15) rate of increase of spacing between successive 
	%		basis functions
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [nU x nB] array where y_ij is the number of spikes at time bin j 
	%			for unit i.
	%		X = [nB x nK] array where X_ijk is the value of covariate k, at time
	%			bin j, for unit i
	%			Note: nK = nU*nK_sp + 2*nK_pos
	%		k = Names of each filter, a [n x 2] cell array in which each row is 
	%			of the form ['filter j name', [idxj1 idxj2 ...]]
	%			Note: The second column lists indices in 1:nK to which the label applies
	%  
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 20;
	%	nK_stm = 6;
	%	a = 10;
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm_network(processed, nK_sp, nK_stm, a);

	if nargin < 4
		a = 15;
	end

	dt_sp = processed.binsize;
	dt_stm = processed.binsize;
	unit = processed.unitidx;

	steps_sp = dt_sp/processed.binsize;
	steps_stm = dt_stm/processed.binsize;

	nB = size(processed.stim,1);
	nU = size(processed.spikes,2);
	nS = size(processed.stim,2);

	T = nK_sp*dt_sp;
	[rcbasis, spbasis, nK_rc] = makeRCBasis(dt_sp, T, a);
	nK = nU*nK_rc + nS*nK_stm;

	data.X = zeros(nB, nK);
	data.k = cell(nU+1,3);
	for i = 1:nU
		data.k{i,1} = ['spike unit' num2str(i)]; 
		data.k{i,2} = (i-1)*nK_rc + (1:nK_rc);
		data.k{i,3} = dt_sp;
	end
	data.k{nU+1,1} = 'stim';
	data.k{nU+1,2} = (1:(nK_stm*nS)) + nU*nK_rc;
	data.k{nU+1,3} = dt_stm;

	%For each unit, add data to X array
	%Make stimulus vector at each timebin
	for j = (nK_sp*steps_sp+1):(nB)
		shist = zeros(nK_rc*nU,1);
		for idx=1:nU 
			%(past) spike history
			shist(((idx-1)*nK_rc+1):(idx*nK_rc)) = project_rc(processed.spikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx), rcbasis);
		end
		%(past) stimulus
		stim = processed.stim(j-nK_stm*steps_stm:steps_stm:j-steps_stm,:);
		stim = reshape(squeeze(stim)', 1, []);
		%Form stim vector
		data.X(j,:) = [shist' stim];
	end
	%Truncate to exclude start and end of recording where spike history 
	%and cursor trajectory aren't well defined
	data.X = data.X((nK_sp*steps_sp+1):end,:); %(nkt+1:end-nkt,:);
	data.y = processed.spikes((nK_sp*steps_sp+1):end,unit)';
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

function [rcbasis, spbasis, nK_rc] = makeRCBasis(dt, T, a)
	%Define basis of raised cosine functions
	if nargin < 3
		a = 15;
	end
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