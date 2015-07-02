function data = filters_sp_stm_network(processed, nK_sp, nK_stm)
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
	%	data = filters_sp_stm_network(processed, nK_sp, nK_stm)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	nK_sp = number of timebins used for spike history filter for all units
	%	nK_stm = number of timebins used for stim filter
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
	%	nK_sp = 6;
	%	nK_stm = 6;
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sp_stm_network(processed, nK_sp, nK_stm);

	dt_sp = processed.binsize;
	dt_stm = processed.binsize;
	unit = processed.unitidx;

	steps_sp = dt_sp/processed.binsize;
	steps_stm = dt_stm/processed.binsize;

	nB = size(processed.stim,1);
	nU = size(processed.spikes,2);
	nS = size(processed.stim,2);

	nK = nU*nK_sp + nS*nK_stm;

	data.X = zeros(nB, nK);
	data.k = cell(nU+1,3);
	for i = 1:nU
		data.k{i,1} = ['spike unit' num2str(i)]; 
		data.k{i,2} = (i-1)*nK_sp + (1:nK_sp);
		data.k{i,3} = dt_sp;
	end
	data.k{nU+1,1} = 'stim';
	data.k{nU+1,2} = (1:(nK_stm*nS)) + nU*nK_sp;
	data.k{nU+1,3} = dt_stm;

	%For each unit, add data to X array
	%Make stimulus vector at each timebin
	for j = (nK_sp*steps_sp+1):(nB)
		shist = zeros(nK_sp*nU,1);
		for idx=1:nU 
			%(past) spike history
			shist(((idx-1)*nK_sp+1):(idx*nK_sp)) = processed.spikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx);
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

end