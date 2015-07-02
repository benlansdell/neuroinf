function model = MLE_glmfit_network(data, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
	%
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%		converged = [nU x 1] array listing 1 if the IRLS converged within iteration limit, 0 if not
	%		conditioned = [nU x 1] array listing 1 if the IRLS did not issue an ill-conditioned warning, 0 if it did
	%
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 6;
	%	nK_stm = 6;
	%	const = 'on';
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm_network(processed, nK_sp, nK_stm);
	%	model = MLE_glmfit_network(data, const);

	if (nargin < 2) const = 'on'; end
	nU = size(data.y,1);
	nK = size(data.X,2);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);
	model.converged = ones(nU,1);
	model.conditioned = ones(nU,1);
	model.nspikes = zeros(nU,1);
	%For each unit, fit a GLM to the torque data
	display(['Fitting GLM by MLE with IRLS. Fitting ' num2str(nU) ' units.'])
	for idx=1:nU 
		[b, dev, stats] = glmfit(data.X,data.y(idx,:),'poisson', 'constant', const);
		%Catch if a warning was raised about badly conditioned matrix
		[warn, warnid] = lastwarn;
		if ~strcmp(warn, '')
	   		switch warnid
        	case 'stats:glmfit:IterationLimit'
        		model.converged(idx) = 0;
        	case 'stats:glmfit:BadScaling'
        		model.conditioned(idx) = 0;
       		end
	    end
	    lastwarn('')
	    %Extract filters fitted...
		model.b_hat(idx,:) = b;	
		model.dev{idx} = dev;
		%Remove residual components since these take up a lot of memory
		model.N = size(stats.resid,1);
		stats = rmfield(stats, {'resid', 'residp', 'residd', 'resida', 'wts'});
		model.stats{idx} = stats;
		model.nspikes = sum(data.y(idx,:));
	end
	%model.logli = ll_network(model, data, 'poisson');
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat]
	end
	display('Done')
